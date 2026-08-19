#!/usr/bin/env wolframscript
(* Unit gate for the codegen EMISSION GUARDS and the stage-handoff contract in
   mathematica/Codegen.m. Symbolic only — no generation, no C++, no compile — so it runs in a plain
   `ctest` in about a second, unlike codegen_regen which takes tens of minutes behind an opt-in flag
   and is therefore not run while editing.

   WHAT IT PINS, and what breaks if it regresses:

   1. $ntCppLeakPatterns names only symbols that EXIST. The list used to hardcode the strings
      "colFacG[", "compileColG[", "compileTInv[", "colourFactorProd[" — and colourFactorProd had
      already been deleted, so the emission backstop was watching for a function that could no
      longer be produced. A rename that updates a definition but not the string list disables the
      guard SILENTLY: the scan still runs, still passes, and catches nothing. Deriving the patterns
      from the symbols makes that impossible, and this test is what proves the derivation is live.

   2. ntExportCpp actually ABORTS on each pattern. A guard that is present but never fires reads as
      coverage while providing none. Checked with $NumTracerVerbose = False, because the whole leak
      scan lives inside an ntLog[] argument and would stop running entirely if ntLog ever acquired
      a Hold attribute — which its own comment used to invite.

   3. ntStageResult refuses a malformed stage record. It is the structural counter to the bug class
      recorded at the head of mkGenerateKernel: a value assigned in an inner scope and read in an
      outer one silently stays an unassigned symbol, and `FreeQ[<unassigned>, Complex]` is
      vacuously True, so "PruneRealTraces" emitted a wrong kernel with no diagnostic.

   Run standalone:  math -noprompt -script numtracer/tests/test_codegen_stages.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

(* 23 fixed checks + one per leak head. Derived, not hardcoded: the head list is meant to GROW as
   new dispatchers appear, and a literal count would turn every such addition into a spurious
   failure — which is how a counted-check guard gets "fixed" by deleting it. Still catches what it
   is for: an Abort escaping mid-file leaves nRan short of this. *)
$nExpected = 23 + Length[$ntLeakHeads];
$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);

(* HoldRest so the expected message is not evaluated to its TEXT before Quiet sees it — see the
   same note in test_ad_param_typing.wls. *)
SetAttributes[abortsQ, HoldRest];
abortsQ[expr_Hold, msg_] := CheckAbort[Quiet[ReleaseHold[expr]; False, msg], True];

Print["=== 1. every head named in $ntCppLeakPatterns is a DEFINED symbol ==="];
(* The desynchronisation this closes had already half-happened: colourFactorProd was in the list and
   deleted from the package. If a rename lands without updating the list, the corresponding name
   here has no DownValues and this fails in ~1 s instead of silently disabling the backstop. *)
Module[{heads = List @@ Map[SymbolName, $ntLeakHeads]},
  check["the head list is non-empty", Length[heads] > 0];
  Do[
    check["head `" <> h <> "` has DownValues",
      DownValues[Evaluate@Symbol["NumTracer`Private`" <> h]] =!= {}],
    {h, heads}];
  check["every head has a matching \"<name>[\" pattern",
    AllTrue[heads, MemberQ[$ntCppLeakPatterns, # <> "["] &]]];

Print["=== 2. ntExportCpp aborts on a leaked head, with verbosity OFF ==="];
(* Verbosity OFF is the point: the scan runs inside an ntLog[] argument, so this is also the
   regression test for "someone gave ntLog a Hold attribute and the backstop stopped running". *)
Module[{tmp = FileNameJoin[{$TemporaryDirectory, "nt_leak_probe.hh"}], saved = $NumTracerVerbose},
  $NumTracerVerbose = False;
  check["a leaked colourFacStr head aborts",
    abortsQ[Hold[ntExportCpp[tmp, "double f(){ return colourFacStr[x]; }"]], MakeNTKernel::cppleak]];
  check["a leaked nt* DSL head aborts",
    abortsQ[Hold[ntExportCpp[tmp, "double f(){ return ntMetric[a,b]; }"]], MakeNTKernel::cppleak]];
  check["a CForm'd Mathematica List aborts",
    abortsQ[Hold[ntExportCpp[tmp, "double f(){ return dressing(List(1,2)); }"]], MakeNTKernel::cppleak]];
  check["a Module/Unique scoped symbol aborts",
    abortsQ[Hold[ntExportCpp[tmp, "double f(){ return tr$2994; }"]], MakeNTKernel::cppleak]];
  check["a private-context symbol aborts",
    abortsQ[Hold[ntExportCpp[tmp, "double f(){ return NumTracer_Private_flavDelta(F1,F2); }"]],
            MakeNTKernel::cppleak]];
  check["Indeterminate aborts",
    abortsQ[Hold[ntExportCpp[tmp, "double f(){ return Indeterminate; }"]], MakeNTKernel::cppleak]];
  (* The complement matters as much: an ordinary identifier containing "nt" must NOT trip it, or
     every kernel using `constant`/`int`/`point` fails to emit. *)
  check["clean C++ containing 'constant'/'int' does NOT abort",
    !abortsQ[Hold[ntExportCpp[tmp, "int f(){ const double constant = 1.0; return (int)constant; }"]],
             MakeNTKernel::cppleak]];
  Quiet @ DeleteFile[tmp];
  $NumTracerVerbose = saved];

Print["=== 3. ntStageResult accepts a well-formed record and refuses a malformed one ==="];
check["a complete record passes through unchanged",
  ntStageResult["probe", {"a", "b"}, <|"a" -> 1, "b" -> 2|>] === <|"a" -> 1, "b" -> 2|>];
check["a MISSING key aborts",
  abortsQ[Hold[ntStageResult["probe", {"a", "b"}, <|"a" -> 1|>]], ntStageResult::keys]];
check["an EXTRA key aborts",
  abortsQ[Hold[ntStageResult["probe", {"a"}, <|"a" -> 1, "b" -> 2|>]], ntStageResult::keys]];
(* THE bug class: a field left as its own unassigned private symbol. It is not Missing and not
   $Failed, it is a Symbol with no value — which is exactly why FreeQ[..., Complex] said True. *)
check["an UNASSIGNED private symbol aborts",
  abortsQ[Hold[ntStageResult["probe", {"a"}, <|"a" -> NumTracer`Private`ntNeverAssignedProbe|>]],
          ntStageResult::unbound]];
check["a Missing value aborts",
  abortsQ[Hold[ntStageResult["probe", {"a"}, <|"a" -> Missing["nope"]|>]], ntStageResult::unbound]];
check["a $Failed value aborts",
  abortsQ[Hold[ntStageResult["probe", {"a"}, <|"a" -> $Failed|>]], ntStageResult::unbound]];
(* An empty list / 0 / False are all legitimate stage outputs and must survive. *)
check["empty and falsy values are legitimate",
  ntStageResult["probe", {"a", "b", "c"}, <|"a" -> {}, "b" -> 0, "c" -> False|>] =!= $Failed];

Print["=== 4. ntkPruneSpec: the PruneRealTraces verdict, and the guards the old bug lacked ==="];
(* THE historic bug: diagData reached this computation as an unassigned symbol, FreeQ[<unassigned>,
   Complex] was vacuously True, every group was judged real, and a kernel shipped with every
   imaginary half dropped. Nothing in the repo could see it — generation was clean, net counts
   unchanged, every oracle green (the pruned kernel computes fewer numbers, all of them right; only
   the DROPPED halves are wrong). These checks are the only place that behaviour is observable.
   groups entries are {diagIdx, ...}; diagData is indexed 1-based by diagIdx+1. *)
check["prune not requested -> empty verdict, whatever the data says",
  ntkPruneSpec[{1, 2}, {{0}, {1}}, True, False, False, True, True]["pruneG"] === {}];
check["a real group is flagged, a complex one is not",
  ntkPruneSpec[{1, 2 + 3 I}, {{0}, {1}}, True, False, True, False, False]["pruneG"] === {True, False}];
(* An unassigned diagData is now impossible (it is a parameter), so the residual failure is a table
   that is too SHORT for the grouping — the same "read the wrong diagram" outcome. Must abort. *)
check["grouping referencing past the end of diagData ABORTS",
  abortsQ[Hold[ntkPruneSpec[{1}, {{0}, {5}}, True, False, True, False, False]],
          mkGenerateKernel::prunedata]];
(* And the bug's exact signature: a COMPLEX flow in which every group came back real. *)
check["complex flow with every group judged real ABORTS",
  abortsQ[Hold[ntkPruneSpec[{1, 2}, {{0}, {1}}, True, False, True, False, False]],
          mkGenerateKernel::pruneall]];
check["the same all-real verdict is FINE when the flow is not complex",
  ntkPruneSpec[{1, 2}, {{0}, {1}}, False, False, True, False, False]["realOnlyG"] === {True, True}];
(* Ordering: when the probe will run, pass 1 must generate UNPRUNED so the probe sees the residual
   it is checking for. Pruning pass 1 is an O(1) kernel error, not a slow path. *)
check["probe will run -> pass 1 emits unpruned",
  ntkPruneSpec[{1, 2 + 3 I}, {{0}, {1}}, True, False, True, True, True]["realOnlyG"] === {}];
check["probe will run is False when offline",
  ntkPruneSpec[{1, 2 + 3 I}, {{0}, {1}}, True, True, True, True, True]["probeWillRun"] === False];

Print[""];
If[nRan =!= $nExpected,
  Print["  FAIL  check count is ", nRan, ", expected ", $nExpected,
        " — a check was added/removed (update $nExpected) or an Abort escaped mid-file"];
  $fail++,
  Print["  PASS  all ", $nExpected, " checks ran (no early Abort)"]];
Print[""];
If[TrueQ[$fail > 0],
  Print["FAILED: ", $fail, " of ", nRan, " check(s)"]; End[]; Exit[1],
  Print["all ", nRan, " checks passed"]; End[]; Exit[0]];
