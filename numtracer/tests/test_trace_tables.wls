#!/usr/bin/env wolframscript
(* Unit gate for the Wolfram side of trace tables ("TraceTables", mathematica/CodegenTables.m):

     1. option normalisation;
     2. scale inference from the frame components: exact homogeneity (degree 1, or 0 for unit
        directions), the external scale unique among the coordinate arguments, ambiguity and failure
        reported, explicit scales checked;
     3. the class traits: masks false exactly on the scales;
     4. the kernel text: every kernel()/post() definition gets a table twin reading ttrN, the stride and
        entries pointer sit after the leading parameters, every kernel() (one per verdict branch) gets a
        tab() with kernel()'s template header and parameters, and nothing else changes.

   Run standalone:  wolfram -script numtracer/tests/test_trace_tables.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);
SetAttributes[abortsQ, HoldRest];
abortsQ[expr_Hold, msg_] := CheckAbort[Quiet[ReleaseHold[expr]; False, msg], True];

Print["=== 1. option ==="];
check["False -> off", ntTraceTableOption[False] === None];
check["Automatic -> default balance, inferred scales, not explicit",
  ntTraceTableOption[Automatic] === <|"Balance" -> $ntTraceTableBalance, "Scales" -> Automatic, "Explicit" -> False|>];
check["True -> every tabulable trace, explicit", ntTraceTableOption[True]["Balance"] === 0. && ntTraceTableOption[True]["Explicit"]];
check["Association -> merged, explicit",
  ntTraceTableOption[<|"Scales" -> {"S0"}|>] === <|"Balance" -> $ntTraceTableBalance, "Scales" -> {"S0"}, "Explicit" -> True|>];
check["unknown key aborts", abortsQ[Hold[ntTraceTableOption[<|"Foo" -> 1|>]], MakeNTKernel::ttopt]];

Print["=== 2. scale inference ==="];
(* a triangle frame: l = l1 (c1, s1 c2, s1 s2, 0), p1 = S0 (x, y, 0, 0) with x, y shape factors of S1 *)
Module[{l1, S0, S1, c1, s1, c2, s2, x, y, usyms, comp, nc},
  usyms = {l1, S0, c1, s1, c2, s2, x, y};
  comp = <|0 -> {l1 c1, l1 s1 c2, l1 s1 s2, 0}, 1 -> {S0 x, S0 y, 0, 0}, 2 -> {l1 c1 - S0 x, l1 s1 c2 - S0 y, l1 s1 s2, 0}|>;
  nc = <|"usyms" -> usyms, "compExpr" -> comp,
         "varFill" -> {"l1", "S0", "cos1", "sqrt(1. - powr<2>(cos1))", "cos2", "sqrt(1. - powr<2>(cos2))", "x(S1)", "y(S1)"}|>;
  With[{s = ntTraceTableScales[nc, "l1", {"S0", "S1", "SPhi"}, Automatic]},
    check["inferred: radial l1 = var 0, external S0 = var 1", s["Radial"] === 0 && s["External"] === 1 && s["ExternalName"] === "S0"]];
  check["explicit {S0} accepted", ntTraceTableScales[nc, "l1", {"S0"}, {"S0"}]["External"] === 1];
  check["explicit {} (no external scale) fails: the momenta need S0",
    FailureQ[ntTraceTableScales[nc, "l1", {"S0"}, {}]]];
  check["explicit non-symbol fails", FailureQ[ntTraceTableScales[nc, "l1", {"S1"}, {"S1"}]]];
  check["radial that is not a symbol of its own fails", FailureQ[ntTraceTableScales[nc, "q", {"S0"}, Automatic]]];
  check["two explicit external scales: not supported", FailureQ[ntTraceTableScales[nc, "l1", {"S0", "S1"}, {"S0", "S1"}]]];
  (* an external vector without a bare scale reads as a unit direction (degree 0); if it is really
     dimensionful, the generator's exact check (table_form) rejects the traces *)
  With[{nc2 = Append[nc, "compExpr" -> <|0 -> {l1 c1, l1 s1 c2, 0, 0}, 1 -> {x, y, 0, 0}|>]},
    check["an external vector without a bare scale is a direction: no external scale",
      ntTraceTableScales[nc2, "l1", {"S0"}, Automatic]["External"] === -1]];
  (* unit directions e = p/|p| ride the frame as degree-0 vectors *)
  With[{nc4 = Append[nc, "compExpr" -> Append[comp, 3 -> {x, y, 0, 0}]]},
    check["unit direction vectors (degree 0) are allowed", ntTraceTableScales[nc4, "l1", {"S0"}, Automatic]["External"] === 1]];
  (* a vector mixing degrees is not homogeneous *)
  With[{nc5 = Append[nc, "compExpr" -> Append[comp, 3 -> {l1 c1 + x, 0, 0, 0}]]},
    check["a vector mixing degree 0 and 1 fails", FailureQ[ntTraceTableScales[nc5, "l1", {"S0"}, Automatic]]]];
  (* a vacuum bubble: only the loop momentum *)
  With[{nc3 = Append[nc, "compExpr" -> <|0 -> {l1 c1, l1 s1, 0, 0}|>]},
    With[{s = ntTraceTableScales[nc3, "l1", {"S0"}, Automatic]},
      check["no external momentum: no external scale", s["External"] === -1 && s["ExternalName"] === None]]];
];

Print["=== 3. traits ==="];
Module[{spec = <|"RadialName" -> "l1", "ExternalName" -> "S0", "BoxNames" -> {"l1", "cos1", "cos2"},
                 "ExtNames" -> {"S0", "S1", "SPhi"}|>, tr},
  tr = ntTraceTableTraits[spec, "DiFfRG::za3_ymf"];
  check["has_tab and tab_size from the traces header",
    MemberQ[tr, "static constexpr bool has_tab = DiFfRG::za3_ymf::tab_size > 0;"] &&
    MemberQ[tr, "static constexpr std::size_t tab_size = DiFfRG::za3_ymf::tab_size;"]];
  check["box mask false on the radial variable only",
    MemberQ[tr, "static constexpr std::array<bool, 3> tab_box_mask{false, true, true};"]];
  check["external mask false on the scale only",
    MemberQ[tr, "static constexpr std::array<bool, 3> tab_ext_mask{false, true, true};"]];
];

Print["=== 4. kernel text ==="];
kernelText = "    template<typename T8, typename T9>
    static KOKKOS_INLINE_FUNCTION auto kernel(const float& l1, const float& cos1, const float& S0, const float& k, const T8& ZA, const T9& Zc)
    {
      float fenv[DiFfRG::ns::nenv];
      DiFfRG::ns::fill(fenv, l1, cos1, S0);
      const auto _a = ZA(l1);
      return _a * DiFfRG::ns::tr0(fenv) + Zc(k) * DiFfRG::ns::tr12(fenv);
    }
    static constexpr std::array<bool, 6> pre_reads{true, false, false, true, true, false};
    template<typename TP_, typename T8, typename T9>
    static KOKKOS_INLINE_FUNCTION auto post(const TP_ &P_, const float& l1, const float& cos1, const float& S0, const float& k, const T8& ZA, const T9& Zc)
    {
      const auto &_a = DiFfRG::device::get<0>(P_);
      float fenv[DiFfRG::ns::nenv];
      DiFfRG::ns::fill(fenv, l1, cos1, S0);
      return _a * DiFfRG::ns::tr0(fenv) + Zc(k) * DiFfRG::ns::tr12(fenv);
    }";
res = ntTableKernelText[kernelText, "DiFfRG::ns", "float", {"l1", "cos1", "S0"}];
kTab = StringCases[res, Shortest["auto kernel_tab(" ~~ __ ~~ "\n    }"]];
pTab = StringCases[res, Shortest["auto post_tab(" ~~ __ ~~ "\n    }"]];
tab = StringCases[res, Shortest["void tab(" ~~ __ ~~ "\n    }"]];
check["the originals are untouched", StringStartsQ[res, kernelText] || StringContainsQ[res, StringTake[kernelText, 300]]];
check["one kernel_tab, one post_tab, one tab", Length[kTab] === 1 && Length[pTab] === 1 && Length[tab] === 1];
check["kernel_tab: entries pointer and stride first",
  StringContainsQ[First[kTab], "auto kernel_tab(const float *row_, const std::size_t stride_, const float& l1,"]];
check["post_tab: after the pre() entry",
  StringContainsQ[First[pTab], "auto post_tab(const TP_ &P_, const float *row_, const std::size_t stride_, const float& l1,"]];
check["every trace reads the table, none the inline trace",
  StringCount[First[kTab], "DiFfRG::ns::ttr0(row_, stride_, fenv)"] === 1 &&
  StringCount[First[kTab], "DiFfRG::ns::ttr12(row_, stride_, fenv)"] === 1 &&
  StringFreeQ[First[kTab], "::tr0(" | "::tr12("] && StringFreeQ[First[pTab], "::tr0(" | "::tr12("]];
check["the twins keep their template headers",
  StringCount[res, "template<typename T8, typename T9>"] === 3 && StringCount[res, "template<typename TP_, typename T8, typename T9>"] === 2];
check["tab(): kernel()'s parameters, fill arguments forwarded",
  StringContainsQ[First[tab], "void tab(double *out_, const std::size_t stride_, const float& l1, const float& cos1, const float& S0, const float& k, const T8& ZA, const T9& Zc)"] &&
  StringContainsQ[First[tab], "DiFfRG::ns::tab_fill(out_, stride_, l1, cos1, S0);"]];
check["pre_reads stays between kernel and post", StringCount[res, "pre_reads"] === 1];
(* a complex kernel defines kernel() once per verdict branch: each branch needs its own tab() *)
With[{k2 = "    #if NT_X_VERDICT == 2\n" <> kernelText <> "\n    #else\n" <> kernelText <> "\n    #endif"},
  With[{r2 = ntTableKernelText[k2, "DiFfRG::ns", "float", {"l1", "cos1", "S0"}]},
    check["verdict branches: one tab() and one kernel_tab() per branch",
      StringCount[r2, "void tab("] === 2 && StringCount[r2, "auto kernel_tab("] === 2 &&
      StringPosition[r2, "void tab("][[1, 1]] < StringPosition[r2, "#else"][[1, 1]] < StringPosition[r2, "void tab("][[2, 1]]]]];

Print[""];
Print[If[$fail === 0, "ALL " <> ToString[nRan] <> " PASSED", ToString[$fail] <> " of " <> ToString[nRan] <> " FAILED"]];
End[];
Exit[If[$fail === 0, 0, 1]];
