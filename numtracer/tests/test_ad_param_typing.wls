#!/usr/bin/env wolframscript
(* Unit gate for the C++ TYPE given to each runtime kernel parameter — ntRuntimeParamType and
   ntAssertADTyped in mathematica/Codegen.m.

   WHY IT EXISTS. An AD-flagged scalar must be emitted `const auto&`, not `const double&`: FunKit's
   MakeParameterString turns "auto" into an abbreviated function template, and that is the only
   reason kernel()/constant() can bind autodiff::real from DiFfRG's integrator_AD twin. The
   ParameterOrder rework (fc7fb4e) silently lost that — a With[] closed one line early, so the
   MemberQ that gates the AD branch tested an UNBOUND symbol, which is False rather than an error,
   and every AD scalar came out `double`.

   WHY IT IS A UNIT TEST AND NOT A FLOW FIXTURE. This particular regression is invisible to every
   other stage that has a test, which is why it reached a user:

     - generation succeeds and prints nothing unusual;
     - net counts are unchanged, so the diagram-count watchdog in regen_check.sh cannot see it;
     - the emitted kernels are NUMERICALLY IDENTICAL, so every compare_* oracle stays green;
     - the ordinary get() compiles and runs correctly;
     - only */src/AD_get.cc — emitted here but compiled by the CONSUMER project — fails.

   On top of that, "ADParams" is reachable only through MakeNTKernelDiFfRG, which no test in this
   repo calls, and not one committed kernel signature under tests/gen contains `const auto&`. So a
   full codegen regeneration is blind to this BY CONSTRUCTION, not by oversight. A symbolic unit
   test is the only place the property is observable without building a DiFfRG consumer.

   WHAT IS PINNED, and what breaks if it regresses:

     1. An AD-flagged scalar types "auto" even though its own ParameterOrder entry says
        "Type" -> "double". This is the exact shape DiFfRG_compat builds (the AD flag lives in a
        separate "AD" key) and the exact shape that regressed. If this loosens, every flow with an
        autodiff twin fails to build downstream.

     2. The AD test runs BEFORE the scalar-name gate. adNames is a subset of scalarNames only as an
        accident of DiFfRG_compat filtering both on Type === "double"; a direct MakeNTKernel caller
        passing ADParams without ScalarParams must still get "auto".

     3. Non-AD parameters are untouched: a plain scalar stays "double", a dressing keeps its own
        entry["Type"], and a bare (non-Association) dressing name falls back to dressTy. This is
        what keeps every committed kernel byte-identical.

     4. Symbol and String spellings of a name agree, on both the parameter side and the ADParams
        side. Hand-written flow files spell names as Symbols; DiFfRG_compat passes Strings.

     5. ntAssertADTyped actually aborts. A guard that silently returns is worse than no guard,
        because it reads as coverage.

   Run standalone:  math -noprompt -script numtracer/tests/test_ad_param_typing.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

(* Total number of check[] calls below. An Abort escaping mid-file would otherwise leave a passing
   tally simply because the remaining checks never executed; asserting the count closes that. *)
$nExpected = 22;
$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);

(* Abort[] is how every guard in Codegen.m reports; CheckAbort turns that into a value so one
   refusal does not tear down the rest of the file. Quiet suppresses only the message we expect.
   HoldRest is load-bearing: without it `MakeNTKernel::adtype` evaluates to the message TEXT before
   Quiet ever sees it, Quiet rejects it (Quiet::mgre) and the abort is not caught — the check then
   reports FAIL against working code. *)
SetAttributes[abortsQ, HoldRest];
abortsQ[expr_Hold, msg_] :=
  CheckAbort[Quiet[ReleaseHold[expr]; False, msg], True];

(* stand-in for the dressTy closure mkGenerateKernel builds from "DressingType"/interpTy. *)
dressTy = Function[nm, "FALLBACK_" <> ntParamName[nm]];

Print["=== 1. ntParamName: Symbol and String spellings collapse to the same key ==="];
check["a String is itself",            ntParamName["d1V"] === "d1V"];
check["a Symbol gives its short name", ntParamName[Global`d1V] === "d1V"];
(* the two spellings MUST agree, or every MemberQ downstream silently misses. *)
check["the two agree",                 ntParamName[Global`d1V] === ntParamName["d1V"]];

Print["=== 2. AD scalars type auto, non-AD scalars stay double ==="];
(* THE regression shape: an Association whose own "Type" is "double", flagged AD only via the
   separate "AD" key that the type field knows nothing about. *)
check["AD scalar, entry Type double -> auto",
  ntRuntimeParamType[<|"Name" -> "d1V", "Type" -> "double", "AD" -> True|>,
    {"d1V", "d2V"}, {"d1V", "d2V", "rhoL"}, dressTy] === "auto"];
check["second AD scalar -> auto",
  ntRuntimeParamType[<|"Name" -> "d2V", "Type" -> "double", "AD" -> True|>,
    {"d1V", "d2V"}, {"d1V", "d2V", "rhoL"}, dressTy] === "auto"];
check["plain scalar -> double",
  ntRuntimeParamType[<|"Name" -> "rhoL", "Type" -> "double"|>,
    {"d1V", "d2V"}, {"d1V", "d2V", "rhoL"}, dressTy] === "double"];
(* empty ADParams is the default (Options[MakeNTKernel]) and must leave everything alone. *)
check["no ADParams at all -> double",
  ntRuntimeParamType[<|"Name" -> "d1V", "Type" -> "double"|>,
    {}, {"d1V"}, dressTy] === "double"];

Print["=== 3. the AD test precedes the scalar-name gate ==="];
(* adNames is a subset of scalarNames only because DiFfRG_compat.m builds both from
   Type === "double". A direct MakeNTKernel caller need not honour that; if the scalar gate ran
   first, this would fall through to entry["Type"] and hand back "double". *)
check["AD name absent from scalarNames -> still auto",
  ntRuntimeParamType[<|"Name" -> "d1V", "Type" -> "double"|>,
    {"d1V"}, {}, dressTy] === "auto"];
check["AD name absent from scalarNames, bare Symbol -> still auto",
  ntRuntimeParamType[Global`d1V, {"d1V"}, {}, dressTy] === "auto"];

Print["=== 4. dressings keep their own type ==="];
check["dressing Association -> its own Type",
  ntRuntimeParamType[<|"Name" -> "ZA", "Type" -> "TexLinearInterpolator1D"|>,
    {"d1V"}, {"d1V"}, dressTy] === "TexLinearInterpolator1D"];
(* the backend-agnostic path (ParameterOrder -> Automatic) maps over BARE names, not Associations,
   and those must reach dressTy — that is how a mixed 1-D/3-D interpolator flow gets per-name
   types. *)
check["bare Symbol dressing -> dressTy",
  ntRuntimeParamType[Global`ZA, {"d1V"}, {"d1V"}, dressTy] === "FALLBACK_ZA"];
check["bare String dressing -> dressTy",
  ntRuntimeParamType["ZA", {"d1V"}, {"d1V"}, dressTy] === "FALLBACK_ZA"];
(* a dressing is never in scalarNames, but if one ever were, the scalar gate must not steal it
   ahead of an AD flag. *)
check["AD flag beats a dressing entry Type",
  ntRuntimeParamType[<|"Name" -> "ZA", "Type" -> "TexLinearInterpolator1D"|>,
    {"ZA"}, {}, dressTy] === "auto"];

Print["=== 5. Symbol/String spellings agree across BOTH sides ==="];
(* four combinations: name spelled as Symbol or String, ADParams entry likewise. All four must
   land on "auto" — a mismatch here is precisely the silent-miss failure mode. *)
check["String name, String adName",
  ntRuntimeParamType[<|"Name" -> "d1V", "Type" -> "double"|>, {"d1V"}, {}, dressTy] === "auto"];
check["Symbol name, String adName",
  ntRuntimeParamType[Global`d1V, {"d1V"}, {}, dressTy] === "auto"];
check["String name, adName normalised from a Symbol",
  ntRuntimeParamType[<|"Name" -> "d1V", "Type" -> "double"|>,
    ntParamName /@ {Global`d1V}, {}, dressTy] === "auto"];
check["Symbol name, adName normalised from a Symbol",
  ntRuntimeParamType[Global`d1V, ntParamName /@ {Global`d1V}, {}, dressTy] === "auto"];

Print["=== 6. ntAssertADTyped ==="];
(* mkParam's shape, reproduced here rather than called, so the guard is tested against the
   Association layout it actually receives. *)
mk[nm_, ty_] := <|"Name" -> nm, "Type" -> ty, "Const" -> True, "Reference" -> True|>;
good = {mk["k", "double"], mk["d1V", "auto"], mk["rhoL", "double"], mk["ZA", "TexLinearInterpolator1D"]};
bad  = {mk["k", "double"], mk["d1V", "double"], mk["rhoL", "double"]};

check["a correctly typed list passes through unchanged",
  ntAssertADTyped[good, {"d1V"}] === good];
check["an empty ADParams list never fires",
  ntAssertADTyped[bad, {}] === bad];
(* an ADParams name that is not a runtime parameter at all must not fire either — MakeNTKernel is
   called with the flow's full ADParams even when a given name did not survive into the signature. *)
check["an ADParams name absent from the list does not fire",
  ntAssertADTyped[good, {"notAParameter"}] === good];
check["an AD name typed double ABORTS",
  abortsQ[Hold[ntAssertADTyped[bad, {"d1V"}]], MakeNTKernel::adtype]];
(* the guard must survive the shape the bug actually produced: SEVERAL AD scalars downgraded at
   once (d1V and d2V together, which is what every FE-potential flow has). *)
check["several AD names typed double ABORTS",
  abortsQ[
    Hold[ntAssertADTyped[{mk["d1V", "double"], mk["d2V", "double"]}, {"d1V", "d2V"}]],
    MakeNTKernel::adtype]];

Print[""];
(* NOT routed through check[]: its argument is evaluated before check[] bumps the counter, so the
   expected number would have to be one less than the total it reports — an off-by-one waiting to be
   "fixed" wrongly by the next person. Asserted directly instead. *)
If[nRan =!= $nExpected,
  Print["  FAIL  check count is ", nRan, ", expected ", $nExpected,
        " — a check was added/removed (update $nExpected) or an Abort escaped mid-file"];
  $fail++,
  Print["  PASS  all ", $nExpected, " checks ran (no early Abort)"]];
Print[""];
(* TrueQ, not a bare comparison: see the note in test_complex_projection.wls — outside this
   context an undefined $fail would leave If unevaluated and the script would exit 0 having
   rendered no verdict at all. *)
If[TrueQ[$fail > 0],
  Print["FAILED: ", $fail, " of ", nRan, " check(s)"]; End[]; Exit[1],
  Print["all ", nRan, " checks passed"]; End[]; Exit[0]];
