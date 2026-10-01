(* Code generation: the imaginary-part probe (ntProbeSource emits it, ntRunProbe builds and runs it
   and writes the verdict header), bounded compiler-log reporting, and diagColPolys, the build-time
   helper that folds group-diagonal dressed colour nets through the C++ engine.
   Loaded by NumTracer.m via ntLoadPart, in the NumTracer`Private` context. *)

(* ---- semantic complexQ: does the imaginary part actually vanish? ---------------------------
   The syntactic `!FreeQ[Diagrams, Complex]` only sees that SOME coefficient carries an i (projector
   i, colour f^abc T^b T^c = (iN/2) T^a). Those i's often pair up across diagrams so the flow is REAL,
   but that depends on the trace VALUES, which are opaque at the Mathematica stage. The probe evaluates
   the generated C++ traces over random frames with smooth real stub dressings and decides.

   The verdict is applied by the PREPROCESSOR: `kernel.hh` carries all three bodies under
   `#if <MACRO> == 2 / #elif == 1 / #else`, and the probe's C++ `main` writes `<MACRO>` into the verdict
   header (`-o <file> -m <MACRO>`). So the probe is an ordinary build step and generation can run
   without a Wolfram kernel.

   "TraceArrayDecl": with CrossTraceCSE the integrand's trace tokens are `tarr[i]` reads, so the probe
   must declare and fill that array exactly as the kernel's coreBlock does. Empty for the per-trace path. *)
Options[ntProbeSource] = {"NPoints" -> 4000, "Tol" -> 1.*^-9, "TraceArrayDecl" -> ""};

(* ---- ntProbeSource: CONTRACT ------------------------------------------------------------------
   Emits the C++ source of the imaginary-part PROBE: a standalone program that evaluates the same
   integrand three ways over @p "NPoints" random kinematic points and prints the residuals the
   verdict is read from.

   In:   integrand   the assembled symbolic integrand (complex),
         args        every kinematic symbol it reads,  fillArgs  the subset the fill needs,
         angleDefs / angleDecls   the named angle temporaries, so the probe TU has them in scope,
         nsHome      the namespace the generated traces live in,
         headerFile  the traces header to include,  drAtoms  the dressing atoms in id order.
   Out:  the probe program text. It is compiled and run by ntRunProbe, which parses its stdout.

   It compares the untouched COMPLEX form with the two real projections (ntPureIntegrand,
   ntRePartIntegrand); the verdict is 2 = Pure, 1 = RePart, 0 = neither (keep it complex).

   Conservative by design: a residual above tolerance resolves to 0 (NaN points are skipped), and
   anything that prevents a verdict (compile, run or parse failure) aborts generation in ntRunProbe.
   A wrong verdict is an O(1) kernel error; an unnecessarily complex kernel is only slower. *)
ntProbeSource[integrand_, args_, fillArgs_, angleDefs_, angleDecls_, nsHome_, headerFile_, drAtoms_List : {}, opts : OptionsPattern[]] :=
  Module[{keepHeads, keepSyms, seedOf, argComb, stub, probeFull, probeProj, probeRePart, probeParams, probePre, fnFull, fnProj, fnRePart, drDecls, drFillArgs, randDecls, callArgs, src, np, tol, distOf},
    np = OptionValue["NPoints"];
    (* A float kernel's roundoff is ~1e-7 relative, so the double tolerance would read noise as a
       surviving Im. 1e-4 is still far below a genuine Im or a dropped term, which are O(1) relative. *)
    tol = If[ntSingleQ[], Max[OptionValue["Tol"], 1.*^-4], OptionValue["Tol"]];
    (* Replace EVERY external real-valued atom (dressing, regulator, named constant) with
       `ntStub(seed_head, hash(args))`, an INDEPENDENTLY-SEEDED pseudo-random real per head. "External":
       not a structural math head, frame arg, known constant, or trace token (raw C++ strings, left
       untouched). Independence is essential: one shared stub would make `dress1[x]-dress2[x]` collapse
       to 0 and mask a real surviving imaginary part. *)
    keepHeads = Alternatives[Plus, Times, Power, Rational, Sqrt, Sin, Cos, Tan, Cot, Sec, Csc, Exp, Log, Abs, Sign, ArcTan, ArcSin, ArcCos, Sinh, Cosh, Tanh, Max, Min, Floor, Ceiling, Mod, Complex, List, Global`ntStub];
    keepSyms = Join[args, First /@ angleDefs, {Pi, E, EulerGamma, Degree, GoldenRatio}];
    (* distinct seed per external head *)
    seedOf[h_] := N[Mod[Hash[SymbolName[h]], 100003] + 7];
    (* 0 for a constant *)
    argComb[a_List] := Total[MapIndexed[#1 * N[GoldenRatio] ^ (#2[[1]] - 1)&, a]];
    stub[e_] := Module[{x},
        x = e //. (h_Symbol)[a___] /; FreeQ[keepHeads, h] && !MemberQ[args, h] :> Global`ntStub[seedOf[h], argComb[{a}]];
        x /. (s_Symbol) /; !MemberQ[keepSyms, s] && FreeQ[keepHeads, s] && Context[s] =!= "System`" && s =!= Global`ntStub :> Global`ntStub[seedOf[s], 0.]
      ];
    probeFull = stub[integrand];
    (* The candidates must be exactly the expressions the Pure and RePart branches EMIT, not proxies:
       e.g. `probeFull /. Complex[a_,b_] :> a` agrees with ntPureIntegrand only for a LINEAR integrand,
       and would certify a multilinear body that drops the −Im(A)·Im(B) legs. Stubbing leaves the trace
       tokens untouched, so the projections see the same tokens they will emit. *)
    probeProj = ntPureIntegrand[probeFull];
    probeRePart = ntRePartIntegrand[probeFull];
    probeParams = (<|"Name" -> SymbolName[#], "Type" -> $ntRealT, "Const" -> True, "Reference" -> True|>)& /@ args;
    (* Dressed kernels: the generated `fill()` takes one extra `double dr_<id>` per dressing atom. The
       probe has no runtime to compute them, so it evaluates each atom with the SAME stubbing as the
       integrand; without them the probe does not compile (too few args to fill). Atoms can reference
       the derived angles, so the angle decls precede them. *)
    drDecls =
      MapIndexed[
        Function[{atom, pos},
          "const " <> $ntRealT <> " dr_" <> ToString[pos[[1]] - 1] <> " = " <> cppFlat[stub[atom]] <> ";"],
        drAtoms];
    drFillArgs = ("dr_" <> ToString[#])& /@ Range[0, Length[drAtoms] - 1];
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
    (* Three more full COEN lowerings of the stubbed integrand; timed separately from the kernel's own
       [prof] body lines because this is pure verdict overhead. *)
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
        (* Generic in the ARGUMENT type: under "ComplexRuntimeProjection" the probe raises a COMPLEX
           denominator (l0 + I muq) to a power. Mirrors numtracer::compute::powr (runtime.hpp). *)
        "template<int N, class T> static inline T powr(T x){ T r=T(1); int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?T(1)/r:r; }\n",
        (* keepHeads KEEPS Max and Min, and unlike the kernel the probe has no `using namespace DiFfRG`
           to supply min/max, so they come from <algorithm> here. *)
        "using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs; using std::min; using std::max;\n",
        With[{c = "std::complex<" <> $ntRealT <> ">"},
          "static inline " <> c <> " fma(const " <> c <> "&a,const " <> c <> "&b,const " <> c <> "&c){return a*b+c;}\n"],
        "template<class T> using complex = std::complex<T>;\n",
        (* Independently-seeded pseudo-random real in [0.4,0.9): same (seed,arg) -> same value (a
           dressing is a function), distinct (seed,arg) -> independent value. Single precision uses a
           SMOOTH stub instead: the hash amplifies its argument ~1e4x, and the compared bodies compute
           that argument in different orders, so float roundoff would give O(1) different dressings
           per body. *)
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
        (* COUNT the RePart points that disagree rather than taking the max. RePart and the complex
           body are associated differently, so an isolated badly-cancelling point is normal roundoff;
           an emitter defect (double-counted or dropped term) is wrong at essentially EVERY point. *)
        "    { double rr = std::fabs(rp-re)/(std::abs(f)+1.0);\n",
        "      if(std::isfinite(rr)){ mrrep=std::max(mrrep, rr); if(rr > " <> ToString[CForm[N[tol]]] <> ") ++nrep; } }\n",
        (* PER-POINT relative measures (mrim, mrdiff): a global max|Im|/max|Re| can let a localized
           Im hide behind a large |Re| at some OTHER point. The +1 floor turns them into an absolute
           test when |Re|, |f| are small. The verdict keys on these; the absolute trio is log-only. *)
        "    if(std::isfinite(im)&&std::isfinite(re)&&std::isfinite(df)){ mim=std::max(mim,std::fabs(im)); mdiff=std::max(mdiff,df); mre=std::max(mre,std::fabs(re));\n",
        "      mrim=std::max(mrim, std::fabs(im)/(std::fabs(re)+1.0)); mrdiff=std::max(mrdiff, df/(std::abs(f)+1.0)); ++ok; } }\n",
        (* Three-way verdict, decided in the C++ so the probe is one self-contained build step:
             0 "Complex"  Im survives                           -> keep the flow complex.
             2 "Pure"     Im=0 AND Complex->Re projection exact -> drop imaginary coefficients.
             1 "RePart"   Im=0 but projection differs           -> a trace is itself complex; take
                                                                   the re/im split.
           No usable points exits nonzero (the caller aborts) rather than baking a guess into a
           committed header. *)
        "  if(ok < 1){ std::fprintf(stderr, \"[probe] no usable points\\n\"); return 2; }\n",
        "  const int verdict = (mrim > " <> ToString[CForm[N[tol]]] <> ") ? 0 : ((mrdiff <= " <> ToString[CForm[N[tol]]] <> ") ? 2 : 1);\n",
        (* the six measures + the RePart outlier count + point count + verdict, on one line. Printed
           BEFORE the RePart gate below so a failure still reports every measure the caller logs. *)
        "  std::printf(\"%.10e %.10e %.10e %.10e %.10e %.10e %ld %ld %d\\n\", mim, mdiff, mre, mrim, mrdiff, mrrep, nrep, ok, verdict);\n",
        (* RePart is an IDENTITY (Re(Σ c·tr) via ntRe/ntIm), so it must reproduce real(probe_full)
           whichever verdict wins; this is the emitter self-check. Gated on the FRACTION of disagreeing
           points: the 1% cut is far from both a badly-conditioned flow (a handful of outliers in 4000)
           and an emitter bug (every point). Errors below the tolerance are invisible by construction. *)
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
   A failing compile can produce a multi-megabyte log, which a Mathematica message renders as an
   opaque `<<N>>` placeholder. So show the head, say how much was dropped, and name the file. Empty
   and missing logs are normal (a run failure often writes nothing to stderr); ReadString returns
   EndOfFile, not "", for an empty file, hence the StringQ guard. *)
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
   header. Any failure ABORTS: the verdict selects a preprocessor branch in a committed header, and
   falling back to the complex body would not bind to the real integrators the DiFfRG scaffold
   declares. Returns the verdict string for the log. *)
ntRunProbe::probefail = "Imaginary-part probe failed: `1`";
ntRunProbe[srcFile_String, tracesDir_String, verdictFile_ : None, macro_ : None] :=
  Module[{cxx = resolveGenCxx[], bin, rc, out, parsed, oflag},
    (* the SOURCE is a committed build input in gen/; the binary and logs are scratch, outside the
       source tree, and unique per invocation so concurrent sessions never share them. *)
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
    (* unique per invocation, so no flow or concurrent session can read another's `.out` *)
    cppFile = FileNameJoin[{$TemporaryDirectory, "ntdiagpoly_" <> StringReplace[CreateUUID[], "-" -> ""] <> ".cpp"}];
    bin = StringReplace[cppFile, ".cpp" -> ""];
    ntExportCpp[cppFile, src];
    (* HEADER_ONLY: sun_value_dressed is a split entry point whose body normally links from
       libNumTracer.a. This helper links nothing, so it needs the bodies inline (cheap: the TU is tiny). *)
    rc = Run[cxx <> " -std=c++20 -O1 -w -DNUMTRACER_HEADER_ONLY=1 -I '" <> includeDir <> "' '" <> cppFile <> "' -o '" <> bin <> "' 2> '" <> bin <> ".cerr'"];
    If[rc =!= 0,
      Print["[diagpoly] compile failed (rc=", rc, "):\n", ntLogHead[bin <> ".cerr"]];
      Abort[]];
    (* Check the run, not just the compile: a crashed helper would otherwise yield res = {} and a
       confusing length mismatch downstream. Delete first so a stale file is never read as output. *)
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
