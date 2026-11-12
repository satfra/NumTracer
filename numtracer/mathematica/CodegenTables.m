(* ::Package:: *)

(* ---- TRACE TABLES ("TraceTables") ----------------------------------------------------------------
   A dressing-free trace is homogeneous in the loop momentum l and the dimensionful external
   coordinate S. Cleared of its projector denominators it is N = Σ_n a_n l^n S^(D-n), whose
   coefficients a_n depend only on the dimensionless kinematics (angles, shape). The generator
   computes them exactly (numtracer/codegen/trace_table.hpp) and emits
     tab_fill(out, stride, <fill args>)   fp64, writes the entries of one table node;
     ttrN(row, stride, f)                 trace N from the node's entries, or a forwarder to trN(f)
                                          where the inline trace is cheaper (the balance);
     tab_size                             entries per node (0: nothing tabulated).
   DiFfRG fills the table once per run over the angle/shape nodes and calls kernel_tab()/post_tab(),
   which this file emits as copies of kernel()/post() reading ttrN instead of trN.

   A dressed or complex trace is tabulated per dressing monomial and real / imaginary part; the
   dressings stay runtime factors. Scope: vacuum flows with a real frame. Anything else falls back
   to the inline kernel (with a log line under Automatic, an error when tables were requested
   explicitly). *)

(* Default balance under Automatic, in flops of inline trace per byte of table entry read. *)
$ntTraceTableBalance = 1.;

MakeNTKernel::ttopt = "\"TraceTables\" -> `1` is not False, True, Automatic or an Association with keys \"Balance\" and/or \"Scales\".";
MakeNTKernel::ttoff = "\"TraceTables\" was requested, but this flow cannot be tabulated: `1`.";

(* -> None (off) | <|"Balance" -> B, "Scales" -> Automatic | {name} | {}, "Explicit" -> bool|> *)
ntTraceTableOption[v_] :=
  Which[
    v === False || v === None, None,
    v === Automatic, <|"Balance" -> $ntTraceTableBalance, "Scales" -> Automatic, "Explicit" -> False|>,
    v === True, <|"Balance" -> 0., "Scales" -> Automatic, "Explicit" -> True|>,
    AssociationQ[v] && SubsetQ[{"Balance", "Scales"}, Keys[v]] &&
      NumericQ[Lookup[v, "Balance", 0.]] && MatchQ[Lookup[v, "Scales", Automatic], Automatic | {___String}],
      <|"Balance" -> N @ Lookup[v, "Balance", $ntTraceTableBalance], "Scales" -> Lookup[v, "Scales", Automatic],
        "Explicit" -> True|>,
    True, Message[MakeNTKernel::ttopt, v]; Abort[]];

(* Is every frame vector homogeneous in the symbols `ids` (0-based Poly var ids), of degree 1 (a momentum)
   or 0 (a unit direction, e.g. e = p/|p| for a projector)? Dimensional consistency then makes every trace
   homogeneous in them. A wrong choice is never silent: table_form re-checks the traces exactly. *)
ntScaleDegreeOneQ[compExpr_Association, usyms_List, ids_List] :=
  AllTrue[Values[compExpr],
    Function[comps,
      With[{degs = DeleteDuplicates[Flatten[(Total[#[[1]][[ids + 1]]]& /@ CoefficientRules[#, usyms])& /@ comps]]},
        Length[degs] <= 1 && SubsetQ[{0, 1}, degs]]]];

(* Does symbol `id` appear in any component? *)
ntSymbolAppearsQ[compExpr_Association, usyms_List, id_Integer] :=
  AnyTrue[Flatten[Values[compExpr]], !FreeQ[#, usyms[[id + 1]]] &];

(* The momentum scales of a frame as Poly var ids: the radial symbol `radialName` and the one external
   coordinate among `candidates` (or none) for which every frame vector is homogeneous of degree 1 or 0;
   an external scale that qualifies is preferred over none (with none, it would be a table axis). `explicit`
   is Automatic (infer) or a list of at most one external name. A scale must be a frame symbol of its
   own (its fill is the bare argument name). Returns <|"Radial", "External", "RadialName",
   "ExternalName"|> (External -1 / None without one), or Failure["TraceTables", <|"Message" -> why|>]. *)
ntTraceTableScales[ncomp_Association, radialName_String, candidates_List, explicit_] :=
  Module[{usyms = ncomp["usyms"], comp = ncomp["compExpr"], idOf, rId, opts, passing},
    idOf[name_] := With[{pos = FirstPosition[ncomp["varFill"], name, None, {1}]},
      If[pos === None, None, pos[[1]] - 1]];
    rId = idOf[radialName];
    If[rId === None,
      Return[Failure["TraceTables", <|"Message" -> "the radial variable " <> radialName <> " is not a frame symbol of its own"|>], Module]];
    opts =
      If[explicit === Automatic,
        Prepend[Select[DeleteCases[idOf /@ candidates, None], ntSymbolAppearsQ[comp, usyms, #]&], None],
        If[Length[explicit] > 1,
          Return[Failure["TraceTables", <|"Message" -> "more than one external scale (" <> StringRiffle[explicit, ", "] <> ") needs bivariate tables, not supported yet"|>], Module]];
        If[explicit === {}, {None},
          With[{id = idOf[First[explicit]]},
            If[id === None,
              Return[Failure["TraceTables", <|"Message" -> "the scale " <> First[explicit] <> " is not a frame symbol of its own"|>], Module]];
            {id}]]];
    passing = Select[opts, ntScaleDegreeOneQ[comp, usyms, DeleteCases[{rId, #}, None]]&];
    If[Length[passing] > 1, passing = DeleteCases[passing, None]];
    Which[
      Length[passing] === 1,
        <|"Radial" -> rId, "External" -> Replace[First[passing], None -> -1], "RadialName" -> radialName,
          "ExternalName" -> If[First[passing] === None, None, ncomp["varFill"][[First[passing] + 1]]]|>,
      passing === {},
        Failure["TraceTables", <|"Message" -> "no choice of external scale among " <> ToString[candidates] <>
          " makes the frame vectors homogeneous together with " <> radialName|>],
      True,
        Failure["TraceTables", <|"Message" -> "the external scale is ambiguous among " <>
          ToString[ncomp["varFill"][[DeleteCases[passing, None] + 1]]] <> "; name it with \"Scales\""|>]]];

(* The per-flow table spec, or None. `args` are the kernel's frame arguments (box variables first);
   `coordNames` the external coordinate argument names. *)
ntTraceTableSpec[o_, ncomp_, args_List, coordNames_, matsubaraSym_] :=
  Module[{tt = o["TraceTables"], cands, why = None, sc},
    If[tt === None, Return[None, Module]];
    cands = Replace[coordNames, Automatic -> {"p"}];
    why =
      Which[
        TrueQ[ncomp["complex"]], "its frame is complex",
        matsubaraSym =!= None, "it is a finite-temperature flow",
        True, None];
    (* DiFfRG reads the external coordinates as the arguments right after the integration variables *)
    If[why === None,
      With[{names = SymbolName /@ args},
        With[{box = TakeWhile[names, !MemberQ[cands, #]&]},
          If[Take[names, UpTo[Length[box] + Count[names, _?(MemberQ[cands, #]&)]]] =!=
               Join[box, Select[names, MemberQ[cands, #]&]],
            why = "its coordinate arguments do not directly follow the integration variables"]]]];
    If[why === None,
      sc = ntTraceTableScales[ncomp, SymbolName[First[args]], cands, tt["Scales"]];
      If[FailureQ[sc], why = sc["Message"]]];
    If[why =!= None,
      If[tt["Explicit"],
        Message[MakeNTKernel::ttoff, why]; Abort[],
        ntLog["[tables] off: ", why]; Return[None, Module]]];
    ntLog["[tables] scales: ", sc["RadialName"], " (radial)",
      If[sc["ExternalName"] === None, "", ", " <> sc["ExternalName"] <> " (external)"], ", balance ", tt["Balance"]];
    Join[sc, <|"Balance" -> tt["Balance"], "Explicit" -> tt["Explicit"],
      "BoxNames" -> TakeWhile[SymbolName /@ args, !MemberQ[cands, #]&],
      "ExtNames" -> Select[SymbolName /@ args, MemberQ[cands, #]&]|>]];

(* ---- kernel.hh: traits and the table entry points --------------------------------------------- *)

(* The class traits: has_tab, tab_size, and which box variables / external coordinates the table
   depends on (all but the scales). `members` as for ntTableKernelText: tab_size is the sum of the
   members' entries. *)
ntTraceTableTraits[spec_, members_List] :=
  With[{mask = Function[names, "{" <> StringRiffle[
          If[# === spec["RadialName"] || # === spec["ExternalName"], "false", "true"]& /@ names, ", "] <> "}"],
        size = StringRiffle[(#["Ns"] <> "::tab_size")& /@ members, " + "]},
    {"static constexpr bool has_tab = " <> size <> " > 0;",
     "static constexpr std::size_t tab_size = " <> size <> ";",
     "static constexpr std::array<bool, " <> ToString[Length[spec["BoxNames"]]] <> "> tab_box_mask" <> mask[spec["BoxNames"]] <> ";",
     "static constexpr std::array<bool, " <> ToString[Length[spec["ExtNames"]]] <> "> tab_ext_mask" <> mask[spec["ExtNames"]] <> ";"}];

(* {start, end} of the member definition whose name token `head` (e.g. "auto kernel(") starts at
   `st`: from its `template<...>` line (past `// clang-format` guards) to its closing brace. *)
ntDefinitionSpan[text_String, chars_List, st_Integer, head_String] :=
  Module[{lineStart = ntLineStart[chars, st], hdr, prev, bOpen},
    hdr = lineStart;
    While[hdr > 1 &&
        StringStartsQ[prev = StringTrim[StringTake[text, {ntLineStart[chars, hdr - 1], hdr - 1}]], "// clang-format"],
      hdr = ntLineStart[chars, hdr - 1]];
    If[hdr > 1 && StringStartsQ[StringTrim[StringTake[text, {ntLineStart[chars, hdr - 1], hdr - 1}]], "template<"],
      lineStart = ntLineStart[chars, hdr - 1]];
    bOpen = First[Select[Range[ntMatchingBracket[chars, st + StringLength[head] - 1, {"(", ")"}] + 1, Length[chars]],
              chars[[#]] === "{" &, 1]];
    {lineStart, ntMatchingBracket[chars, bOpen, {"{", "}"}]}];

(* The entries of member `m` within a table row: the row pointer, offset by the entries of the members
   before it. *)
ntTableRowOf[m_, base_String] :=
  If[m["Offset"] === "0", base, base <> " + (" <> m["Offset"] <> ") * stride_"];

(* A definition of kernel() or post() turned into its table form: the entries pointer and stride after
   the leading parameters, and every trace call of a tabulated member reading that member's entries. *)
ntTableDefinition[fn_String, members_List, realT_String] :=
  Module[{entry = "const " <> realT <> " *row_, const std::size_t stride_, ", out},
    out = StringReplace[fn, {
        "auto kernel(" -> "auto kernel_tab(" <> entry,
        "auto post(const TP_ &P_, " -> "auto post_tab(const TP_ &P_, " <> entry}, 1];
    Fold[
      Function[{t, m},
        StringReplace[t,
          RegularExpression[StringReplace[m["Ns"], "::" -> "\\:\\:"] <> "\\:\\:tr(\\d+)\\(fenv\\)"] :>
            m["Ns"] <> "::ttr$1(" <> ntTableRowOf[m, "row_"] <> ", stride_, fenv)"]],
      out, members]];

(* tab(out, stride, <kernel parameters>) from a kernel() definition: the same template header and
   parameters, filling each tabulated member's entries of one table node from the fill arguments. *)
ntTableFillFn[fn_String, members_List, fillArgNames_List] :=
  Module[{kpos = First[First[StringPosition[fn, "auto kernel(", 1]]], chars = Characters[fn], open, close, prefix,
          params},
    open = kpos + StringLength["auto kernel"];
    close = ntMatchingBracket[chars, open, {"(", ")"}];
    prefix = StringTake[fn, kpos - 1];
    params = StringTake[fn, {open + 1, close - 1}];
    prefix <> "void tab(double *out_, const std::size_t stride_, " <> params <> ")\n    {\n" <>
      StringJoin[("      " <> #["Ns"] <> "::tab_fill(" <> ntTableRowOf[#, "out_"] <> ", stride_, " <>
        StringRiffle[fillArgNames, ", "] <> ");\n")& /@ members] <> "    }"];

(* The kernel text with the kernel_tab()/post_tab() copies inserted after the originals, and tab()
   after each kernel() (one per preprocessor branch of a complex kernel's verdict). `members` lists the
   tabulated trace namespaces in table order, each <|"Ns" -> namespace, "Offset" -> C++ expression for
   the entries before it|> ("0" for the first); trace calls of other namespaces stay inline. *)
ntTableKernelText[text_String, members_List, realT_String, fillArgNames_List] :=
  Module[{chars = Characters[text], spans, out = text, shift = 0},
    spans =
      SortBy[
        Join[
          ({#, "auto kernel("}& /@ (First /@ StringPosition[text, "auto kernel("])),
          ({#, "auto post("}& /@ (First /@ StringPosition[text, "auto post("]))],
        First];
    Do[
      Module[{span = ntDefinitionSpan[text, chars, s[[1]], s[[2]]], fn, add},
        fn = StringTake[text, span];
        add = "\n" <> ntTableDefinition[fn, members, realT];
        If[s[[2]] === "auto kernel(", add = add <> "\n" <> ntTableFillFn[fn, members, fillArgNames]];
        out = StringInsert[out, add, span[[2]] + 1 + shift];
        shift += StringLength[add]],
      {s, spans}];
    out];

(* ---- kernel groups ---------------------------------------------------------------------------- *)

(* The table spec of a kernel group from its members' specs (`specs`, None for a member without a
   table) and trace namespaces `nss`, or None. The tabulated members must agree in their scales and
   table axes; the result carries them in table order as "Members" (see ntTableKernelText). *)
ntGroupTableSpec[specs_List, nss_List, groupName_String, explicit_] :=
  Module[{tab = Select[Range[Length[specs]], specs[[#]] =!= None &], keys = {"RadialName", "ExternalName", "BoxNames", "ExtNames"}},
    If[tab === {}, Return[None, Module]];
    If[!SameQ @@ (KeyTake[specs[[#]], keys]& /@ tab),
      If[TrueQ[explicit],
        Message[MakeNTKernelGroupDiFfRG::differ, groupName, "their trace-table scales"]; Abort[],
        ntLog["[tables] group ", groupName, " off: the members differ in their trace-table scales"];
        Return[None, Module]]];
    Append[specs[[First[tab]]],
      "Members" -> MapIndexed[
        <|"Ns" -> nss[[#1]], "Offset" ->
          If[First[#2] === 1, "0", StringRiffle[(nss[[#]] <> "::tab_size")& /@ Take[tab, First[#2] - 1], " + "]]|>&,
        tab]]];
