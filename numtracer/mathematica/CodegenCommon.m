(* ::Package:: *)
(* CodegenCommon.m — shared helpers of the code generator: verbose logging, C++ literal/table text,
   the single emission chokepoint ntExportCpp (leak scan), the stage hand-off contract, MakeNTKernel
   messages, and net-builder chunking. Loaded by NumTracer.m via ntLoadPart, in NumTracer`Private`.

   Code generation is split over Codegen*.m, loaded in this order: Common (helpers), Nets (Lorentz /
   colour / Dirac factor emission), Frames (momentum components and fill symbols), RealProjection
   (real/imaginary projection of the integrand), Build (compiler, include/lib, manifest), Probe
   (imaginary-part probe), Generator (the build-time generator program), Kernel (the kernel header).

   NumTracer owns only the tensor part: each component is contracted numerically to a polynomial
   (MPoly) and Horner-lowered to straight-line C++. Everything scalar (coefficients, CSE, boilerplate)
   goes through FunKit's COEN emitter. The seam: each trace result is a C++ identifier that appears as
   a string placeholder in the integrand, which CppForm emits verbatim (CExpression[a_String] := a).

   Files of one generation: the generator program (gen_<ns>.cpp, its gen_<ns>_u<k>.cpp units, the
   _nets.hh declarations, the _pch.hh header), the traces header it prints when run
   (<Name>_kernels.hh), and the kernel header the consumer includes (<Name>_kernel.hh). *)

(* Verbose diagnostics ([prof], [cse], [probe], [diagpoly], [time]) go through ntLog and are silent
   unless $NumTracerVerbose (env NT_GEN_VERBOSE=1); errors and "wrote:"/"unchanged:" always print.
   Delayed (:=) so SetEnvironment works after the package is loaded.

   ntLog is HoldAll: its arguments evaluate ONLY when verbose, so it must never wrap load-bearing
   work (a guard placed inside would silently not run). Bind the work outside, log only the timing:
       With[{ntT = First @ AbsoluteTiming[ <the work> ]}, ntLog["[prof] ...: ", ntT, " s"]];
   tests/test_codegen_stages.wls asserts ntExportCpp's leak abort fires with verbosity off. *)
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
   Dense flows carry hundreds of millions of sub-term scalars; mixed exact/precision-17 numbers do not
   pack (~125 vs 16 bytes each). cppNum prints a machine double exactly, so values survive, but the
   per-net merge now rounds in double: literals can differ from the precision-17 path by a few ulp
   (byte-identical when every scalar is exactly representable). NT_GEN_EXACT_SCALARS=1 restores the
   precision-17 path, the one the committed reference kernels come from. *)
$ntExactScalars := Environment["NT_GEN_EXACT_SCALARS"] === "1";
ntPackCx[l_List] :=
  If[$ntExactScalars || !VectorQ[l, NumberQ],
    l,
    Developer`ToPackedArray[N[l] + 0. I]];

(* ---- integer-table text ----------------------------------------------------------------------
   The emitted generator is dominated by flat integer tables; listable IntegerString over the packed
   array is ~2x faster than ToString /@ list.
   TRAP: IntegerString DROPS THE SIGN (IntegerString[-5] is "5"), which would silently miswrite an
   index. Hence the Min >= 0 gate and the ToString fallback. *)

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
   `harvest[]` gives the distinct values in id order. That is deliberately the SAME order as
   `DeleteDuplicates` over the flattened column, so carrying key columns as integer ids keeps the
   emitted bytes unchanged. `Internal`Bag` because appending to a plain list is O(n^2). *)

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
   EVERY generated file goes through here; nothing else may Export generated source.
   An expression that reaches the emitter un-lowered gets ToString'd verbatim into the file: at best
   a cryptic compile error far from the cause, at worst text that is valid C++ and compiles into a
   silently wrong kernel. This textual scan is the backstop for that whole class; the structural
   guards upstream (tleak/colleak/eagernn) give better messages and stay. *)

(* Every private compiler/dispatcher whose UNEVALUATED form could be ToString'd into generated C++.
   Held, so naming a symbol here cannot evaluate it. Derived from the symbols, not written as strings,
   so a rename cannot silently disable the scan; tests/test_codegen_stages.wls asserts every head here
   still has DownValues. *)
$ntLeakHeads = Hold[
    colourFacStr, compileColour, compileColourSum, compileLorentz, compileLorentzBody,
    splitColourGroups, compileDirac, chunkLorentz, lorentzNetStr, lorentzElemStr,
    diracSlotStr, dressedSlotStr];

$ntCppLeakPatterns =
  Join[
    List @@ Map[Function[ntLeakSym, SymbolName[Unevaluated[ntLeakSym]] <> "[", HoldFirst], $ntLeakHeads],
  {
    (* any un-lowered nt* DSL head; the word boundary keeps `constant`, `int`, `point` from matching *)
    RegularExpression["(?<![A-Za-z0-9_])nt[A-Z][A-Za-z0-9]*\\["],
    "Indeterminate",
    "DirectedInfinity",
    "ComplexInfinity",
    "$Failed",
    "Missing[",
    (* A CForm'd Mathematica List: the generic signature of an unresolved consumer-side head (a
       FunKit `dressing[...]`, anything the flow forgot to map), which cannot be enumerated. Lowered
       C++ builds aggregates with braces and has no identifier `List`. *)
    RegularExpression["(?<![A-Za-z0-9_])List\\("],
    (* A scoped symbol (`rho$1767`, `tr$2994`, `ntRad$3`): nullary, so the List( rule misses it, and
       GCC/Clang ACCEPT `$` in identifiers, so it can compile into a wrong kernel. Producers include
       ntSplitRealImag's imaginary-unit stand-in, ntProjectIntegrand's `Unique["tr$"]` placeholders
       and TensorBases' TBUnique dummy indices. No legitimate emitted identifier contains `$`. *)
    RegularExpression["(?<![A-Za-z0-9_$])[A-Za-z][A-Za-z0-9]*\\$[0-9]+"],
    (* A symbol from any PRIVATE context: CForm prints NumTracer`Private`flavDelta[F1,F2] as
       `NumTracer_Private_flavDelta(F1, F2)` — no List( tail, no $nnn, and COEN hoists it like a
       legitimate dressing lookup. A private symbol is never part of an emitted interface. Matches
       `_Private_`, not `NumTracer_`, because the user-supplied kernel "Name" may carry a prefix. *)
    RegularExpression["(?<![A-Za-z0-9_])[A-Za-z][A-Za-z0-9]*_Private_"]}];

(* The scan is ONE StringPosition sweep with the pattern list as alternatives (2.6x faster than one
   sweep per pattern; StringContainsQ with the list is not faster). The pattern that hit is then
   recovered from the 300-character context window the message quotes. *)
ntExportCpp[file_, text_] := (
    (* the scan runs unconditionally: it is bound in the With, and only its timing goes into the
       HoldAll ntLog *)
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
    (* OpenWrite does not create parent directories, and a fresh checkout may lack flows/<name>/ or gen/ *)
    Module[{dir = DirectoryName[file]},
      If[StringQ[dir] && dir =!= "" && !DirectoryQ[dir],
        CreateDirectory[dir, CreateIntermediateDirectories -> True]]];
    (* WriteString, not Export[..., "Text"]: same bytes, without Export's re-encoding cost *)
    Module[{st = OpenWrite[file, CharacterEncoding -> "UTF8"]},
      WriteString[st, text];
      Close[st]]);

(* ---- stage hand-off contract ---------------------------------------------------------------
   A value declared in an outer Module but assigned inside an inner one that also declares it stays
   an unassigned Symbol; downstream reads are then silently wrong (`FreeQ[<unassigned>, Complex]` is
   vacuously True). Every generation stage returns an Association through here, so a field it forgot
   to assign fails at the hand-off, naming stage and field. Empty lists, 0 and False pass; only
   "never assigned" (or Missing/$Failed) is refused. *)

ntStageResult::keys = "Stage `1` returned the keys `2`, but its contract declares `3`. A stage must return exactly the fields it promises — a missing one means a code path forgot to assign it (the failure this guard exists for), an extra one means the contract is stale.";

ntStageResult::unbound = "Stage `1` returned field `2` as `3`, which is not a value (an unassigned private symbol, Missing, or $Failed). An unassigned symbol makes every downstream read silently wrong. Usually the field was assigned inside an inner Module/With that also declared it.";

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

MakeNTKernel::cppleak = "`1`: the generated source still contains un-lowered Mathematica — the text `2` appears in it. Writing it would produce a file that either fails to compile or, worse, compiles into a silently wrong kernel. This means some expression reached the emitter without being turned into C++; the fragment above should identify which. Offending file:\n`3`";

(* max chars of net-builder elements packed into ONE emitted function (see ntChunkDef). A single
   braced-init of a dense net is one huge basic block (slow compile, OOM) and materialises every
   temporary in one full-expression (runtime stack overflow); size-bounded helpers avoid both and let
   the unit bin-packer spread them across TUs. Kept well under $ntUnitChars so chunks stay packable. *)

$ntDefChunk = 60000;

(* target chars per generator TU, and the hard cap on how many TUs a flow may split into. The cap must
   stay above what the size target asks for, or units silently grow back past the target. *)

$ntUnitChars = 250000;

$ntUnitCap = 512;

(* Split a list of C++ element strings into consecutive runs of ~$ntDefChunk characters, counting `sep`
   separator characters per element. Cumulative chars / chunk size is nondecreasing, so equal keys form
   contiguous runs. *)
ntSplitByChars[xs_List, sep_Integer] :=
  SplitBy[Transpose[{xs, Ceiling[Accumulate[(StringLength /@ xs) + sep] / $ntDefChunk]}], Last][[All, All, 1]];

(* Emit ONE net builder `ret name()`. Small builders are a single braced-init. Large ones become
   `name_c<k>(o)` helpers that append in order, plus an assembler `name()` calling them, so the
   vector is element-for-element identical. Two dedup paths come first, because these tables are
   often redundant row-wise: all rows equal -> the fill constructor; few distinct rows -> a distinct
   table plus an index run, taken only when it shrinks the source.
   Returns {defs, decls}: the defs are ordinary top-level defs the bin-packer may scatter across
   units, the decls go into the shared header where the cross-unit calls resolve. *)
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
          (* the distinct table name_u() is called from the assembler, which may land in another
             unit, so it gets its own forward declaration below (ntChunkDef only declares the
             helpers it generates, never its own entry point) *)
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
          (* A flat int table is emitted as `static const int a[] = {...}` + one insert per chunk,
             not one `o.push_back(...); ` per element (~3x less source, far faster -O0 compile).
             Tested on the distinct elements `u`, which is exact and cheap. *)
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
   An oversized main() breaks compilers per FUNCTION: endless -O1 compiles, a runtime stack frame
   beyond 8 MB, and cc1plus stack overflows on deeply nested initializers. Bounded top-level chunk
   functions avoid all three; values and order are unchanged.
   Returns {definitionsString, callExpression}. *)

ntBigTableFns[name_String, ret_String, elems_List] :=
  If[elems === {},
    {"static " <> ret <> " " <> name <> "(){ return {}; }\n", name <> "()"},
    Module[{chunks, nChunks},
      chunks = ntSplitByChars[elems, 2];
      nChunks = Length[chunks];
      (* each chunk is a flat braced-init in its own function (push_back per element costs ~14
         chars each); a single chunk skips the concatenation wrapper *)
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
