(* ::Package:: *)
(* Code generation. NumTracer owns ONLY the tensor part — contracting each
   component numerically to a polynomial (MPoly) that is Horner-lowered to a real
   straight-line kernel. Everything scalar (the dressing/regulator coefficients, CSE,
   powr<n>, the function/class/header boilerplate, clang-format, write-if-changed)
   is delegated to FunKit's mature COEN emitter (CppForm / MakeCppFunction /
   MakeCppClass / MakeCppHeader / WriteCodeToFile), which produces the flat
   straight-line kernel form.

   The seam: each component's scalar result is bound to a C++ identifier in the
   kernel-body PREAMBLE we emit; that identifier appears as a *string placeholder*
   in the Mathematica integrand `Σ coeff_i × trace_i`, which FunKit's CppForm emits
   verbatim (CExpression[a_String] := a) while CSE-ing the coefficients around it. *)
(* Generation goes through the numeric matrix-product backend (`lorentzNetStr`/`compileLorentz` + the
   generator below); generated kernels are validated against FormTracer (FORM) oracles in the
   test suite. *)
(* Files of one generation: the build-time generator program (gen_<ns>.cpp, its gen_<ns>_u<k>.cpp
   units, the _nets.hh declarations and the _pch.hh header), the straight-line traces header it
   prints when run (<Name>_kernels.hh), and the kernel header the consumer includes (<Name>_kernel.hh),
   which fills the fundamental symbols and calls the generated traces. *)
(* Verbose-diagnostics gate. The profiling / CSE / probe / timing traces below ([prof], [cse],
   [probe], [diagpoly], [time]) are emitted through ntLog and stay SILENT unless this flag is set —
   so a normal generation run is quiet. Genuine "[NumTracer] ERROR" aborts and the "wrote:"/
   "unchanged:" file messages are always printed (plain Print). To see the diagnostics, set
   NumTracer`Private`$NumTracerVerbose = True before generating.

   HoldAll: ntLog evaluates its arguments ONLY when verbose, so it must NEVER wrap a side-effecting
   computation. Bind such work in a With and pass only the timing in:

       With[{ntT = First @ AbsoluteTiming[ <the work> ]}, ntLog["[prof] ...: ", ntT, " s"]];

   This is not a style preference. The attribute was added only after every call site was audited,
   because the file used to do the opposite — the ntExportCpp LEAK SCAN and DSL.m's checkLabels
   guard both ran inside an ntLog argument, and they worked only because ntLog was an ordinary
   function. Adding HoldAll while either was still there would have deleted a correctness guard from
   every non-verbose run and left a green test suite. tests/test_codegen_stages.wls pins the
   ntExportCpp half of that (it asserts the abort still fires with verbosity OFF). *)
(* Default silent, but env-controllable so a headless/CI run can turn the [cse]/[prof]/[time]
   diagnostics on without editing a .wls — mirrors the C++ side's NT_GEN_PROFILE. The density guard in
   tests/gen/regen_check.sh needs the [cse] sub-terms line, which is emitted through ntLog.
   Delayed (:=), like every other env-derived flag here, so `SetEnvironment` works from a .wls that
   has already loaded the package — with an immediate `=` the value latched at Get[] time and the
   documented recipe silently did nothing. *)
$NumTracerVerbose := ntEnvFlag["NT_GEN_VERBOSE"];

SetAttributes[ntLog, HoldAll];
ntLog[args___] := If[TrueQ[$NumTracerVerbose],
    Print[args]];

(* WolframKernel resident set (MB) from /proc — the RSS the OOM killer sees, not MemoryInUse[]
   (allocated). Returns 0 where /proc is absent. *)
ntWolframRssMB[] := Quiet @ Check[
    Module[{s = ReadString["/proc/self/statm"]},
      If[StringQ[s], ToExpression[StringSplit[s][[2]]] * 4096. / 1048576., 0.]],
    0.];

(* a C++ double literal at full precision (exact rationals -> doubles). *)

cppNum[x_] := ToString[CForm[If[MachineNumberQ[x], SetPrecision[x, 17], N[x, 17]]]];

(* ---- sub-term scalars as PACKED machine complex ------------------------------------------------
   The expansion below produces one scalar per sub-term, and the dense finite-T flows have tens to
   hundreds of millions of them (a finite-T/finite-mu four-quark lambda4L2: 396 M). They used to be
   carried as the precision-17 numbers that N[num, 17] (slot options) and the branch scalars
   produce, mixed with exact integers, so no per-net list could pack: measured ~125 bytes per
   scalar against 16 packed. Together with packing the two slotCombo id columns (24 -> 8 bytes),
   a sub-term costs ~56 bytes instead of ~190 (measured on that flow: 23 GB after the expansion,
   where the unpacked expansion had passed 36 GB at 85% of the nets).

   Packing them as machine doubles keeps the emitted VALUES: cppNum prints a machine number as
   SetPrecision[x, 17], the exact 17-significant-digit decimal of that double, which a C++ compiler
   parses back to the identical double. What changes is where rounding happens: the per-net merge
   sums shared sub-terms in double instead of at 17 digits (a few ulp), and the emitted literals are
   no longer byte-identical to the precision-17 path. Where every scalar is exactly representable
   (e.g. ZA/ZAPre of a finite-T QCD tree) the output IS byte-identical. NT_GEN_EXACT_SCALARS=1
   restores the precision-17 path exactly -- the one the committed reference kernels come from. *)
$ntExactScalars := Environment["NT_GEN_EXACT_SCALARS"] === "1";
ntPackCx[l_List] :=
  If[$ntExactScalars || !VectorQ[l, NumberQ],
    l,
    Developer`ToPackedArray[N[l] + 0. I]];

(* ---- integer-table text ----------------------------------------------------------------------
   The emitted generator is dominated by flat integer tables: 99.9% of ZAAqbq2's 25.7 MB main TU is
   ntSidxU/ntDscU/sdrU literals (6.0 M numeric tokens). `ToString /@ list` walks them one downvalue
   at a time; `IntegerString` is listable and runs over the packed array in one kernel call —
   measured 8.2 s -> 4.1 s on 6.3 M integers, output identical.

   THE TRAP: IntegerString DROPS THE SIGN. IntegerString[-5] is "5", not "-5". Every table here is
   an index or an exponent and so non-negative today, but a table that ever grew a negative entry
   would be silently miswritten — an off-by-a-sign index into a trace table, i.e. exactly the
   compiles-fine-but-wrong failure this file's guards exist for. Hence the Min >= 0 gate (a packed
   Min is C-speed) and the ToString fallback; the guard costs nothing measurable. *)

ntIntStrs[l_List] :=
  If[l === {},
    {},
    If[VectorQ[l, IntegerQ] && Min[l] >= 0,
      IntegerString[l],
      ToString /@ l]];

(* one braced C++ initialiser row from a list of integers *)

ntIntRow[l_List] := "{" <> StringRiffle[ntIntStrs[l], ","] <> "}";

(* ---- first-appearance interner --------------------------------------------------------------
   Returns {intern, harvest}: `intern[v]` gives v's 0-based id, minting a new one on first sight;
   `harvest[]` gives the distinct values in id order. This is deliberately the SAME order a
   `DeleteDuplicates` over the flattened column would produce, which is what lets the emitter carry
   its key columns as integers from the start without moving a single emitted byte — see the
   emitNumericGenerator note on integer key columns.
   `Internal`Bag` for the value list: appending to a plain list is O(n^2) at the scale this runs. *)

ntMkIntern[] :=
  Module[{idx = <||>, bag = Internal`Bag[], n = 0},
    (* Lookup evaluates its default only on a miss *)
    {Function[v, Lookup[idx, Key[v], idx[v] = n; Internal`StuffBag[bag, v]; n++]],
     Function[Null, Internal`BagPart[bag, All]]}];

(* CppForm, flattened to a single line — the fill formulas are emitted INSIDE C++ string
   literals, so CppForm's line-wrapping (newlines) must be collapsed or they break the string. *)

cppFlat[e_] := StringReplace[FunKit`CppForm[e], {"\n" -> " ", "\r" -> " ", "\t" -> " "}];

(* The real type of EMITTED code: "double", or "float" under "ComputeType" -> "float". MakeNTKernel
   Blocks it (together with FunKit's literal precision) for one generation; every emitter string that
   names a real type or a real literal reads it, so a float kernel carries no double that would
   promote its arithmetic back. The generator itself always computes in double. *)
$ntRealT = "double";

ntSingleQ[] := $ntRealT === "float";

ntZeroLit[] := If[ntSingleQ[], "0.f", "0.0"];

(* ---- the single emission chokepoint --------------------------------------------------------
   EVERY generated file goes through here. Nothing else may call Export on generated source.

   The failure this exists for: an expression that reached the emitter WITHOUT being lowered to C++
   gets ToString'd verbatim into the file. It has happened three times — the ZAAqbq metric leak
   (`colourFacStr[ntMetric[...], <|...|>]`), a degenerate Gram emitting `return Indeterminate;`, and the
   four-quark Fierz flavour sum (`colourFacStr[Plus[...], <|...|>]`). Each time the symptom appeared
   layers away from the cause: a clang syntax error naming a column of a 7000-character line, or —
   the dangerous variant — a leak whose text happens to BE valid C++ and compiles into a silently
   wrong kernel.

   A textual assertion catches the whole class at the one place it must pass through, regardless of
   which upstream dispatcher grew a hole. Structural guards upstream (tleak/colleak/eagernn) give
   better messages and should stay; this is the backstop that cannot be forgotten.

   `nt` heads are matched with a word boundary so ordinary C++ identifiers containing "nt"
   (`constant`, `int`, `point`) do not trip it. *)

(* Every PRIVATE compiler/dispatcher whose UNEVALUATED form would be ToString'd into generated C++.
   Held, not a plain list, so naming a symbol here cannot evaluate it (several are defined further
   down this file, and one is a dispatcher whose bare name would match its own fallthrough rule).

   Derived from the SYMBOLS rather than written out as strings, because the string spelling is what
   rots: the list used to carry a literal "colourFactorProd[" for a function that had already been
   deleted, so the backstop was watching for something that could no longer be produced — and a
   rename that missed a string would have disabled it just as silently. tests/test_codegen_stages.wls
   asserts every head here still has DownValues, so a rename that forgets this list fails in ~1 s. *)
$ntLeakHeads = Hold[
    colourFacStr, compileColour, compileColourSum, compileLorentz, compileLorentzBody,
    splitColourGroups, compileDirac, chunkLorentz, lorentzNetStr, lorentzElemStr,
    diracSlotStr, dressedSlotStr];

$ntCppLeakPatterns =
  Join[
    List @@ Map[Function[ntLeakSym, SymbolName[Unevaluated[ntLeakSym]] <> "[", HoldFirst], $ntLeakHeads],
  {
    RegularExpression["(?<![A-Za-z0-9_])nt[A-Z][A-Za-z0-9]*\\["],
    (* any nt* DSL head, un-lowered *)
    "Indeterminate",
    "DirectedInfinity",
    "ComplexInfinity",
    "$Failed",
    "Missing[",
(* A CForm'd Mathematica List. This is the generic signature of a CONSUMER-side symbolic head that
   no rule ever resolved — the head itself is arbitrary (a FunKit `dressing[...]`, a `GammaN[...]`,
   anything the flow file forgot to map), so it cannot be enumerated, but any such head printed by
   CForm carries its argument lists as `List(...)`. Catching it here turns a leak that otherwise
   costs a full emission plus a failed compile into an immediate, named failure.
   `List(` cannot arise from legitimately lowered code: the emitted C++ builds every aggregate with
   braces, and no runtime/support identifier is spelled `List`. *)
    RegularExpression["(?<![A-Za-z0-9_])List\\("],
(* A Mathematica SCOPED symbol — `Module`/`Block` auto-renaming (`rho$1767`) or `Unique` (`tr$2994`,
   `ntRad$3`, `ffslash$12`). These are NULLARY, so they carry no `List(...)` argument tail and the
   generic rule above cannot see them; and, unlike a stray head, they print as a bare identifier
   that GCC and Clang both ACCEPT ($ in an identifier is a documented extension). So the failure is
   an "undefined identifier" only as long as the name happens not to collide with something
   declared — otherwise it compiles into a silently wrong kernel.
   Three known producers: `ntSplitRealImag`'s local stand-in for the imaginary unit (see there),
   `ntProjectIntegrand`'s multilinear `Unique["tr$"]` placeholders, and
   TensorBases' `TBUnique` dummy indices escaping through the closure returned by
   TB3PToS0S1SPhi / TB3PToS0as (Kinematics.m) when they land in a tensor slot no rule rewrites.
   No legitimate emitted identifier contains `$`: verified zero `$` characters across all 630
   committed generated .hh/.cpp under tests/gen. *)
    RegularExpression["(?<![A-Za-z0-9_$])[A-Za-z][A-Za-z0-9]*\\$[0-9]+"],
(* A symbol from a PRIVATE Mathematica context. CForm renders the context marks as underscores, so
   NumTracer`Private`flavDelta[F1,F2] prints as `NumTracer_Private_flavDelta(F1, F2)`. Neither
   generic rule above can see it: it carries no `List(...)` tail (its arguments are bare index
   symbols) and no `$nnn` (FunKit names an internal index with Unique["F"], which yields `F45`,
   not `F$45` — Routing.m:497). And COEN hoists it into `const auto _interpN = ...`, indistinguishable
   from a legitimate dressing lookup, so nothing upstream objects either. That is how the
   uncontracted fundamental-flavour delta reached a kernel.
   A private symbol is BY CONSTRUCTION not part of any emitted interface, so its appearance in
   generated text is always a leak, whichever package it escaped from — this also covers
   FunKit`Private`* and TensorBases`Private`*. Matching `_Private_` rather than `NumTracer_` is
   deliberate: the kernel "Name" option is user-supplied, so a package prefix could be a legitimate
   identifier. Verified zero matches across all 661 committed generated .hh/.cpp under tests/gen. *)
    RegularExpression["(?<![A-Za-z0-9_])[A-Za-z][A-Za-z0-9]*_Private_"]}];

(* ONE pass, not thirteen. The scan runs over EVERY generated file, and on a dense flow that is tens
   of megabytes (ZAAqbq2: 45.9 MB across 83 files) — pattern-at-a-time meant 13 full sweeps, five of
   them with the regex engine. StringPosition takes the pattern LIST as alternatives and matches them
   in a single sweep: measured 37.7 s -> 14.4 s on a 309 MB string (2.6x). `StringContainsQ` with the
   same list is NOT the faster form (measured 39 s) — do not "simplify" this to it.
   Which pattern hit is recovered afterwards from the 300-character context window, which is small
   enough that re-testing all thirteen there is free. The context is what the message quotes anyway:
   the matched token alone rarely identifies which structure it came from, and the emitted lines are
   thousands of characters wide. *)
ntExportCpp[file_, text_] := (
(* THE SCAN RUNS UNCONDITIONALLY. It is bound here and only its TIMING is logged — the scan must not
   sit inside the ntLog[] call. ntLog is HoldAll, so its arguments evaluate only when verbose: a scan
   placed there would silently not run in every non-verbose generation. Keep load-bearing work
   outside ntLog. *)
    With[{ntT =
        First @
          AbsoluteTiming[
            Module[{pos = StringPosition[text, $ntCppLeakPatterns, 1], ctx, culprit},
              If[pos =!= {},
                ctx = StringTake[text, {Max[1, pos[[1, 1]] - 150], Min[StringLength[text], pos[[1, 2]] + 150]}];
                culprit = FirstCase[$ntCppLeakPatterns, q_ /; StringContainsQ[ctx, q] :> q, "<unidentified>"];
                Message[MakeNTKernel::cppleak, "ntExportCpp", StringTake[text, First[pos]], file <> "\n  pattern: " <> ToString[culprit] <> "\n  context: ..." <> ctx <> "..."];
                Abort[]]]]},
      ntLog["[prof] ntExportCpp leak scan (", Round[StringLength[text] / 1048576.], " MB, ",
        FileNameTake[file], "): ", ntT, " s"]];
    (* ensure the target directory exists — a fresh checkout may have neither flows/<name>/ nor gen/
       yet, and OpenWrite does not create parents (it fails instead). Cheap and idempotent. *)
    Module[{dir = DirectoryName[file]},
      If[StringQ[dir] && dir =!= "" && !DirectoryQ[dir],
        CreateDirectory[dir, CreateIntermediateDirectories -> True]]];
(* WriteString, not Export[...,"Text"]: byte-identical output (verified) but Export re-encodes the
   whole string through its converter stack, which is measurable on a 26 MB main TU. *)
    Module[{st = OpenWrite[file, CharacterEncoding -> "UTF8"]},
      WriteString[st, text];
      Close[st]]);

(* ---- stage hand-off contract ---------------------------------------------------------------
   The structural counter to the bug class recorded at the head of mkGenerateKernel: a value
   DECLARED in an outer Module but ASSIGNED inside an inner one silently stays an unassigned symbol
   — not Missing, not $Failed, just a Symbol with no value. Nothing complains, and the reads that
   follow are quietly wrong: `FreeQ[<unassigned>, Complex]` is vacuously True, which is how
   "PruneRealTraces" -> True once emitted a kernel with every group wrongly pruned.

   Every extracted generation stage returns its outputs as an Association and routes it through here,
   so a field a stage forgot to assign fails AT THE HAND-OFF, naming the stage and the field, instead
   of surfacing hundreds of lines later as a wrong number. Cheap enough to run unconditionally: a
   handful of key comparisons per stage, against generation runs measured in seconds to minutes.

   Empty lists, 0 and False are legitimate stage outputs and pass; only "never assigned" is refused. *)

ntStageResult::keys = "Stage `1` returned the keys `2`, but its contract declares `3`. A stage must return exactly the fields it promises — a missing one means a code path forgot to assign it (the failure this guard exists for), an extra one means the contract is stale.";

ntStageResult::unbound = "Stage `1` returned field `2` as `3`, which is not a value: an unassigned private symbol, Missing, or $Failed. An unassigned symbol is the dangerous case — it is not an error in Mathematica, so every downstream read of it silently produces a wrong answer (FreeQ[<unassigned>, Complex] is vacuously True; that is how PruneRealTraces emitted a wrong kernel). Almost always the field was assigned inside an inner Module/With that also DECLARED it, so the assignment never reached this scope.";

ntStageResult[label_String, keys_List, a_Association] :=
  Module[{bad},
    If[Sort[Keys[a]] =!= Sort[keys],
      Message[ntStageResult::keys, label, Keys[a], keys]; Abort[]];
    bad = Select[keys,
      With[{v = a[#]},
        MissingQ[v] || v === $Failed ||
          (Head[v] === Symbol && Context[v] === "NumTracer`Private`" && ! ValueQ[v])] &];
    If[bad =!= {},
      Message[ntStageResult::unbound, label, First[bad], a[First[bad]]]; Abort[]];
    a];

MakeNTKernel::eagernn = "compileLorentz: an eagerly-summed structure has a NON-NUMERIC per-structure scalar coefficient, which the emitted add(...) cannot carry (each summand is scaled by a numeric literal). The sum should have been distributed by expandBridges (DSL.m) or collected into an ntDressedNum. Offending sum:\n`1`";

MakeNTKernel::tleak = "compileLorentz: un-lowered TENSOR structure reached the scalar fallthrough — it would be CForm'd into C++ as a bare scalar with its indices silently dropped (this is what caused the ZAAqbq metric leak). Every tensor head must be handled by lorentzNetStr or one of the Power/Times/Plus branches. Offending structure:\n`1`";

MakeNTKernel::colleak = "compileColour: a factor of a CONSTANT SU(N) component matched none of the colourFacStr head rules, so it would be ToString'd into the generated C++ as raw Mathematica (e.g. `colourFacStr[Plus[...], <|...|>]`), which does not compile. Every factor of a constant component must be one of the six group heads (ntSUNf/ntSUNDeltaAdj/ntSUNT/ntSUNDeltaFund/ntSUNDiag{Fund,Adj}). A Plus here means a colour/flavour sum reached the flat compiler instead of compileColourSum; any other head means the \"Constant\" head list in DSL.m analyseDiagram has drifted. Offending factor:\n`1`";

(* NB: no backquoted code fragments in this string — a backquoted word is a StringForm SLOT, so
   quoting an identifier that way makes the message itself fail to format (StringForm::sfr). *)

MakeNTKernel::colrest = "splitColourGroups: after splitting a branch into its colour product and its Lorentz/Dirac remainder, SU(N) head(s) are still present in the REMAINDER. The remainder goes to the Lorentz-only lorentzNetStr/compileDirac, which has no rule for a group head, so it would be CForm'd into the generated C++ as raw Mathematica (lorentzNetStr[ntSUNDeltaFund[...], <|...|>]). The split is Cases/DeleteCases at LEVEL 1, so it only sees group heads that are BARE factors — one buried inside a Power (a CLOSED colour loop: deltaFund[N,i,j]^2) or a surviving Plus slips through. Offending head(s):\n`1`\nRemainder:\n`2`";

MakeNTKernel::colpow = "compileColour: a colour/flavour SUM raised to the integer power `1`. Expanding it by repetition would duplicate the summands' index labels, so the same label would appear on 2k tensors and the et/SUNNet contraction would silently mis-pair them into a wrong number — and checkLabels has already run by this point, so nothing downstream would notice. (This is the colour analogue of NumTrace::bridgepow.) Refusing instead. Offending base:\n`2`";

MakeNTKernel::tokleak = "ntProjectIntegrand: a scoped Mathematica symbol survived the real/imaginary projection of the integrand and would be CForm'd into the kernel as a bare C++ identifier (which both GCC and Clang ACCEPT, so it can compile into a silently wrong kernel rather than failing). Two producers: the local stand-in for the imaginary unit in ntSplitRealImag, which survives whenever Coefficient could not treat the coefficient as a polynomial in it (ntIiSafeQ is meant to have refused that case first), and — on the multilinear route — a trace-token placeholder that was not substituted back out. Either way the token-degree routing or the nesting guard missed a summand, most likely an integrand shape that is neither a Plus of Times nor covered by ntSplitTokenPart. Offending symbol(s):\n`1`";

(* NB: no backquoted code fragments in these two strings either — see the note above. *)

MakeNTKernel::cplxnest = "ntProjectIntegrand: the integrand contains an explicit Complex NESTED below a head the real/imaginary split cannot traverse (typically a denominator, i.e. Power with a negative exponent — the finite-density shape l0 + I muq). The split substitutes the imaginary unit by a real stand-in and reads off its coefficients, which is faithful ONLY where the coefficient is a polynomial in that stand-in. For anything else Coefficient power-series expands and returns the leading term, so 1/(a + I b) projects to 1/a: the imaginary part of the denominator is DELETED, silently, and the kernel is a wrong number that compiles. Non-rational heads (Sqrt, Exp) are worse still — the whole expression comes back as the stand-in's zeroth-order part, leaking the local symbol into the emitted C++. Refusing instead. `1` offending subexpression(s):\n`2`";

MakeNTKernel::toknest = "ntRePartLinear: a trace token in this summand is not a BARE factor of it — it sits inside a Power, or below some other head. The token-degree routing classified the summand as linear, but neither the factor-level scan here nor the Coefficient extraction it replaced can see such a token, so the whole summand would be dropped from the projected integrand: a missing term, with nothing downstream to notice. Offending summand:\n`1`";

MakeNTKernel::cppleak = "`1`: the generated source still contains un-lowered Mathematica — the text `2` appears in it. Writing it would produce a file that either fails to compile or, worse, compiles into a silently wrong kernel. This means some expression reached the emitter without being turned into C++; the fragment above should identify which. Offending file:\n`3`";

MakeNTKernel::adtype = "ntRuntimeParamType: `1` runtime parameter(s) named in ADParams were NOT typed auto, so they are emitted as const double& instead of const auto&. That is a SILENT defect here and a failure 20 minutes away in the consumer: an auto parameter is what makes the emitted function an abbreviated template, and it is the only reason kernel()/constant() can bind autodiff::real. Typed double, the flow generates cleanly, the net counts are unchanged, the kernels are numerically identical and the ordinary get() compiles and runs — only the AD twin AD_get.cc fails to instantiate, in a project this repo never builds. That is exactly how the ParameterOrder rework shipped the regression. Cause is almost always that the AD name did not survive the ToString/SymbolName normalisation, or that ADParams carries a name that is not a runtime parameter at all. Offending name(s) and the type each got:\n`2`";

(* max chars of net-builder elements packed into ONE emitted function (see ntChunkDef below). A net
   builder is normally emitted as a single braced-init list of its sub-term elements:
     std::vector<std::vector<DChainTok>> dch<i>(){ return {E1, ..., En}; }
   On a dense flow that one expression grows to megabytes / tens of thousands of elements. At -O0 the
   back end handles it as ONE basic block with thousands of live temporaries, so a single such TU can
   cost minutes and several GB of RSS — the direct cause of both a very long compile AND an OOM, since
   several such units run in parallel. Worse, the braced-init materialises every element temporary in
   one full-expression, so a big enough builder overflows the default stack at RUNTIME. Chunking the
   element list into size-bounded helper functions fixes all three, and — since each helper is its own
   top-level def — lets the unit bin-packer spread them across TUs, so no single unit carries a whole
   giant net. Kept well under $ntUnitChars so chunks stay packable. *)

$ntDefChunk = 60000;

(* target chars per generator TU, and the hard cap on how many TUs a flow may split into. The cap must
   stay above what the size target asks for, or units silently grow back past the target. *)

$ntUnitChars = 250000;

$ntUnitCap = 512;

(* Emit ONE net builder, chunking it when its element list is too big to sit in a single function.
   Small builders keep the old single-def form byte-for-byte. Big ones become

     void dch<i>_c0(std::vector<std::vector<DChainTok>>& o){ o.push_back(E1); o.push_back(E2); ... }
     ...
     std::vector<std::vector<DChainTok>> dch<i>(){
       std::vector<std::vector<DChainTok>> o; o.reserve(n); dch<i>_c0(o); ...; return o; }

   Order is preserved, so the assembled vector is element-for-element identical to the braced-init
   one — and strictly less work at runtime, since a braced-init-list copies every element while
   push_back moves the temporary. Returns {defs, decls}: the helpers go into `allDefs` as ordinary
   top-level defs (so the bin-packer may scatter them across units), and their forward declarations
   go into the shared decl header, which is where the cross-unit calls resolve. *)

(* Two dedup fast paths precede the generic chunking. They exist because these tables are hugely
   redundant ROW-WISE, which the per-element chunker cannot see (measured 2026-08-08 on the emitted
   generator sources):

     table     rows    distinct    share of source (za3_147 / za4)
     sdn0      33618          1        26% / 4%     -> every row is literally `DiracNet{}`
     sln0      33618         19        22% / 70%
     sdslR0    33618      33601        35% / 4%     -> genuinely distinct, correctly left alone

   Both paths reproduce the vector element-for-element; only the SOURCE TEXT shrinks. `sdslR0` is
   the control that keeps this honest: the guard below rejects it on its own merits, because a
   33601/33618 table would pay for an index vector and save nothing. That guard is why no
   opt-out is needed — the dedup is applied only where it demonstrably pays. *)

(* Split a list of C++ element strings into consecutive runs of ~$ntDefChunk characters, counting `sep`
   separator characters per element. Cumulative chars / chunk size is nondecreasing, so equal keys form
   contiguous runs. *)
ntSplitByChars[xs_List, sep_Integer] :=
  SplitBy[Transpose[{xs, Ceiling[Accumulate[(StringLength /@ xs) + sep] / $ntDefChunk]}], Last][[All, All, 1]];

ntChunkDef[name_String, ret_String, elems_List] :=
  Module[{u, tot},
    u = DeleteDuplicates[elems];
    tot = Total[StringLength /@ elems];
    Which[
      Length[elems] <= 1 || tot <= $ntDefChunk,
        {{ret <> " " <> name <> "(){ return {" <> StringRiffle[elems, ", "] <> "}; }"}, ""},

      (* --- every row identical: the fill ctor, and the payload is emitted ONCE --- *)
      Length[u] === 1,
        {{ret <> " " <> name <> "(){ return " <> ret <> "(" <> ToString[Length[elems]] <> ", " <> First[u] <> "); }"}, ""},

      (* --- few distinct rows: distinct table + an index run --- *)
      (* Conservative: compares payloads only, charging the index to the new form while giving the
         old form no credit for the ~14 chars/element of `o.push_back(...); ` it also pays. *)
      (* one index entry costs its digits plus a comma *)
      Length[u] < Length[elems] && Total[StringLength /@ u] + (StringLength[ToString[Length[u]]] + 1) * Length[elems] < tot,
        Module[{uDefs, uDecl, pos, xs, xchunks, xnc, xdefs, xdecl, n},
          pos = Lookup[AssociationThread[u -> Range[Length[u]] - 1], elems];
          n = ToString[Length[elems]];
          {uDefs, uDecl} = ntChunkDef[name <> "_u", ret, u];
          (* the distinct table is reached from the assembler below, which the bin-packer may put in
             another unit, so it needs its own forward declaration (ntChunkDef only declares the
             _c helpers it generates, never its own entry point). *)
          xs = ntIntStrs[pos];
          xchunks = ntSplitByChars[xs, 1];
          xnc = Length[xchunks];
          xdefs = MapIndexed["void " <> name <> "_x" <> ToString[#2[[1]] - 1] <> "(std::vector<int>& o){ static const int a[] = {" <> StringRiffle[#1, ","] <> "}; o.insert(o.end(), a, a + " <> ToString[Length[#1]] <> "); }"&, xchunks];
          xdecl = StringJoin[Table["void " <> name <> "_x" <> ToString[k - 1] <> "(std::vector<int>&);\n", {k, xnc}]];
          {
            Join[uDefs, xdefs,
              {ret <> " " <> name <> "(){ " <> ret <> " u = " <> name <> "_u(); std::vector<int> x; x.reserve(" <> n <> "); " <>
                 StringJoin[Table[name <> "_x" <> ToString[k - 1] <> "(x); ", {k, xnc}]] <>
                 ret <> " o; o.reserve(" <> n <> "); for (int i : x) o.push_back(u[i]); return o; }"}],
            uDecl <> ret <> " " <> name <> "_u();\n" <> xdecl}],

      True,
        Module[
          {intElems, chunks, nChunks, defs},
(* Is this a flat table of integer literals? If so the generic path below emits `static const int
   a[] = {...}` + one insert instead of one `o.push_back(...); ` per element — the SAME dense form
   the index run further down already uses. The boilerplate it drops is ~14 characters per element,
   which on a dressed flow is not a rounding error: ZAAqbq2's `sdchR0` alone was 7.6 MB of the 20 MB
   of net-builder units. Values, order and the resulting vector are unchanged; measured on a 200k
   table, 4.18 MB -> 1.38 MB of source and a -O0 compile of 14.2 s -> 0.44 s.
   Tested on the DISTINCT elements (`u`), not on `elems`: these tables have millions of entries drawn
   from a few hundred distinct indices, so the check is free where it matters and still exact. The
   `ret` gate keeps it off every non-int table (DiracNet/NetVal/DSlotOpt/...). *)
          intElems = ret === "std::vector<int>" && AllTrue[u, StringMatchQ[#, ("-" | "") ~~ DigitCharacter ..]&];
          chunks = ntSplitByChars[elems, 2];
          nChunks = Length[chunks];
          defs =
            If[intElems,
              MapIndexed["void " <> name <> "_c" <> ToString[#2[[1]] - 1] <> "(" <> ret <> "& o){ static const int a[] = {" <> StringRiffle[#1, ","] <> "}; o.insert(o.end(), a, a + " <> ToString[Length[#1]] <> "); }"&, chunks],
              MapIndexed["void " <> name <> "_c" <> ToString[#2[[1]] - 1] <> "(" <> ret <> "& o){ " <> StringJoin[("o.push_back(" <> # <> "); ")& /@ #1] <> "}"&, chunks]];
          {Append[defs, ret <> " " <> name <> "(){ " <> ret <> " o; o.reserve(" <> ToString[Length[elems]] <> "); " <> StringJoin[Table[name <> "_c" <> ToString[k - 1] <> "(o); ", {k, nChunks}]] <> "return o; }"], StringJoin[Table["void " <> name <> "_c" <> ToString[k - 1] <> "(" <> ret <> "&);\n", {k, nChunks}]]}
        ]
    ]];

(* ---- big literal tables as TOP-LEVEL functions, not braced-inits inside main() ---------------
   One oversized main() is a compiler hazard in three separate ways, all measured on a dense flow
   (hPhiL, four-quark basis, 9.94 MB main TU): clang++ -O1 did not finish in an hour; the -O0
   fallback then SEGFAULTED on entry, its stack frame exceeding the 8 MB default; and g++ -O1 died
   with "internal compiler error: Segmentation fault" after 21 s at 1.79 GB — a cc1plus stack
   overflow on the deeply-nested initializers, which `ulimit -s unlimited` alone works around.
   All three are per-FUNCTION effects, so chunking the tables into ordinary top-level functions
   fixes them at once, and the emitted values are unchanged (push_back in the same order, which also
   moves each temporary instead of copying it out of a braced-init).
   Returns {definitionsString, callExpression}. *)

ntBigTableFns[name_String, ret_String, elems_List] :=
  If[elems === {},
    {"static " <> ret <> " " <> name <> "(){ return {}; }\n", name <> "()"},
    Module[{chunks, nChunks},
      chunks = ntSplitByChars[elems, 2];
      nChunks = Length[chunks];
(* Each chunk is a FLAT braced-init returned by its own function, not one push_back per element.
   The push_back form costs ~14 characters of boilerplate per element, which on the flat index
   vectors (tens of thousands of bare integers) inflated the TU by ~47% and cost 67% more compile
   time — the braced-init keeps the original text density while still bounding each function.
   A single chunk skips the concatenation wrapper entirely, so small tables are emitted exactly as
   they were before. *)
      {
        If[nChunks === 1,
          "static " <> ret <> " " <> name <> "(){ return {" <> StringRiffle[First[chunks], ","] <> "}; }\n",
          StringJoin[
            MapIndexed["static " <> ret <> " " <> name <> "_c" <> ToString[#2[[1]] - 1] <> "(){ return {" <> StringRiffle[#1, ","] <> "}; }\n"&, chunks],
            "static " <> ret <> " " <> name <> "(){ " <> ret <> " o; o.reserve(" <> ToString[Length[elems]] <> ");\n" <>
              StringJoin[Table["  { " <> ret <> " c = " <> name <> "_c" <> ToString[k - 1] <> "(); o.insert(o.end(), std::make_move_iterator(c.begin()), std::make_move_iterator(c.end())); }\n", {k, nChunks}]] <>
              "  return o; }\n"]],
        name <> "()"}]];

(* the same over every net of one family (dnet/lnet/dch/dsl); returns {flatDefs, declString}. *)

ntChunkDefs[prefix_String, ret_String, elemLists_List] := If[elemLists === {},
    {{}, ""},
    Module[{r = MapIndexed[ntChunkDef[prefix <> ToString[#2[[1]] - 1], ret, #1]&, elemLists]},
      {Flatten[r[[All, 1]]], StringJoin[r[[All, 2]]]}]];
