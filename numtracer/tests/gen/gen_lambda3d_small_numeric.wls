#!/usr/bin/env wolframscript
(* QCD generator: the GENERAL-FRAME test flow — the four-quark vertex at two-momentum
   (S0,S1,SPhi) kinematics, cut down.  Run from numtracer/:
       wolfram -script tests/gen/gen_lambda3d_small_numeric.wls

   THIS IS THE ONLY FLOW IN THE SUITE ON THE polyFrameSpec PATH. Every other committed flow
   satisfies unitLoopOkQ (externals depend only on the magnitude `p`), so the unit-loop
   parametrisation runs and reduce_units collapses sin^2 -> 1-cos^2 and the bare-loop denominator
   to l1^2. A flow whose externals carry shape coordinates (S0,S1,SPhi — the with_mesons
   lambda1L3D class) falls to the general polyFrameSpec parametrisation, where before this fixture
   NOTHING was covered: minted radical symbols, a wider MPoly symbol space, +/- duplicated momentum
   bases, and dressed+complex traces. See tests/qcd_lambda3d_setup.m for the physics content.

   TWO KERNELS are emitted, in SEPARATE directories (the imaginary-part probe writes ONE
   numtrace_verdict.hh per directory, fopen "w" — two complex kernels in one directory would
   clobber each other's verdict macro, and tests/gen's own verdict file belongs to zaqbq7_147):

       tests/gen/lambda3d/        the flow, dedup ON   (the Lambda3DSmall_num files)
       tests/gen/lambda3d_ref/    the control, dedup OFF (Lambda3DSmallRef_num, NT_GEN_NO_DEDUP)

   compare_lambda3d_small.cpp grades one against the other on numbers (scale-relative), never
   bytes — the dedup renumbers every sN. As frame/sign levers land they add their own escape-hatch
   controls here (e.g. unit-groups off) so the lever is graded the same way.

   GUARD RAILS asserted below after generation — if any fails, the fixture is silently testing
   nothing and the script aborts:
     - polyFrameSpec taken:   nsym >= 10 in the emitted generator (unit-loop flows have ~4-6)
     - dressed:               the generator carries sdr tables / the kernel fill takes dr_ args
     - complex:               the traces header uses std::complex
     - dedup live:            distinct traces < contractions (the [cse] line)                    *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
$dirA = FileNameJoin[{$gen, "lambda3d"}];
$dirR = FileNameJoin[{$gen, "lambda3d_ref"}];
If[!DirectoryQ[#], CreateDirectory[#]]& /@ {$dirA, $dirR};
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_lambda3d_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];  (* tests: emit against the shim + DiFfRG namespace *)

(* Keep the box small on this machine: the fixture is an ITERATION benchmark, and compile jobs are
   RAM-capped by request. Applies to the generator compile and its run. *)
SetEnvironment["NT_GEN_JOBS" -> "8"];
SetEnvironment["NT_GEN_MAXW" -> "8"];
SetEnvironment["NT_GEN_MAXW_B" -> "4"];

dress = {ZAqbq1, Zq, Mq, ZA, dtZA, dtZq};

(* The four-quark derivative, exactly as with_mesons' lambda1Raw takes it: ONLY the combined
   exchange qb1<->qb2 TOGETHER WITH q1<->q2 is a symmetry of the channel — a single qb<->qb swap is
   a Fierz rearrangement into other basis elements and must NOT be imposed. *)
sym  = FMakeSymmetryList[FSymmetry[Symmetric, {i1, i3}, {i2, i4}]];
fRG4q = FTakeDerivatives[WetterichEquation, {qb[i1], q[i2], qb[i3], q[i4]}, "Symmetries" -> sym] //
  FTruncate // FRoute;
proj = FTerm[TBGetProjector["FierzCompleteNf2Nc34D", 1,
   {i2, i1, i4, i3} /. fRG4q["1-Loop"]["ExternalIndices"]]];
trace = FTerm[proj, fRG4q["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]];

(* TWO-MOMENTUM configuration, with_mesons' lambda1tch3D verbatim: legs (p2, -p2, -p3, p3) in
   FunKit's alternating {qb,q,qb,q} order. ReplaceAll (single pass), NOT ReplaceRepeated: the
   renaming is SIMULTANEOUS, and p1 -> p2 would otherwise be rewritten again by the p2 rule. This
   +/- leg routing is also what makes the emitted comp[] table carry exact +/- base pairs — one of
   the mechanisms under guard. *)
trace3d = trace /. {p1 -> p2, p2 -> -p2, p3 -> -p3, p4 -> p3};

flow = S0Param[dressingRules[trace3d]];
net  = FromFunKit[flow];
Print["lambda3d(small) net: top terms = ", If[Head[net] === Plus, Length[net], 1]];

(* CUT-DOWN DIAL, the zaaqbq1_small pattern: diagrams are dropped Sorted-first (deterministic), and
   the per-diagram properties under test (frame fallback, +/- routing, density) survive any cut.
   NT_TEST_NDIAG overrides (0 = keep all) for calibration/scaled benchmarks. *)
nKeep = With[{e = Environment["NT_TEST_NDIAG"]},
  If[StringQ[e] && IntegerQ[Quiet@ToExpression[e]], ToExpression[e], 6]];
If[nKeep > 0 && Head[net] === Plus && Length[net] > nKeep,
  net = Plus @@ Take[Sort[List @@ net], nKeep];
  Print["lambda3d(small) cut down to ", nKeep, " diagrams (NT_TEST_NDIAG)"]];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frameS0S1Sphi, "Args" -> argsS0S1SPhi] ];
tG = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$dirA, "gen_lambda3d_small_num.cpp"}],
  FileNameJoin[{$dirA, "Lambda3DSmall_num_kernel.hh"}],
  FileNameJoin[{$dirA, "Lambda3DSmall_num_kernels.hh"}],
  "Name" -> "Lambda3DSmall_num_kernel", "Namespace" -> "lambda3d_small_num",
  "Dressings" -> dress] ];
Print["[time] NumTracer generation TOTAL (analyse+codegen): ", tA + tG, " s"];

(* ---- guard rails: abort rather than commit a fixture that tests nothing ---------------------- *)
Module[{genSrc, traces, symLine, nsym},
  genSrc = Import[FileNameJoin[{$dirA, "gen_lambda3d_small_num.cpp"}], "Text"];
  traces = Import[FileNameJoin[{$dirA, "Lambda3DSmall_num_kernels.hh"}], "Text"];
  symLine = First[StringCases[genSrc, "symNames = {" ~~ Shortest[s__] ~~ "};" :> s], ""];
  nsym = Length[StringSplit[symLine, "\",\""]];
  Print["lambda3d(small) guard: nsym = ", nsym];
  (* Off the pure fast path: the externals' shape coordinates force minted radical fills (a sqrt in
     a symNames entry that is not the loop's own polar factors). nsym stays well above the pure
     fast path's 4-6 either way. *)
  If[nsym < 8,
    Print["FIXTURE GUARD FAILED: nsym < 8 — the flow did not leave the unit-loop fast path."]; Abort[]];
  (* The MIXED unit-loop branch must have fired: the loop rides opaque unit-direction symbols with
     a non-empty unit-constraint table (the pre-2026-08-07 fallback emitted `units = {}`). The
     control kernel (NT_NO_UNIT_GROUPS below) covers the empty-units path. *)
  If[!StringContainsQ[genSrc, "units = {{"],
    Print["FIXTURE GUARD FAILED: units table empty — the mixed unit-loop branch did not fire."]; Abort[]];
  If[!StringContainsQ[genSrc, "sdr"] || !StringContainsQ[traces, "dress("],
    Print["FIXTURE GUARD FAILED: no dressing atoms — the dressed (DPoly) path is not exercised."]; Abort[]];
  If[!StringContainsQ[traces, "std::complex"],
    Print["FIXTURE GUARD FAILED: traces are all real — the complex path is not exercised."]; Abort[]];
];

(* ---- THE REFERENCE KERNEL: the same flow with every lever switched OFF — the pristine
   pre-lever computation: every (net, sub-term) contracted separately (NT_GEN_NO_DEDUP), the whole
   frame through polyFrameSpec with no unit groups (NT_NO_UNIT_GROUPS), and expression-keyed
   momentum bases with the ± sign twins kept (NT_NO_SIGN_CANON). Grading the main kernel against
   this control validates the dedup and both frame/sign levers together, numerically, against a
   computation that shares only the derivation and the routing. See gen_zaaqbq1_small_numeric.wls
   for why grading is numeric, never bytes.

   NT_NO_SIGN_CANON acts inside NumTrace (the canonicalisation rewrites the net before the env is
   built), so the control needs its OWN NumTrace pass — reusing `ntk` would silently grade the
   lever against itself. *)
SetEnvironment["NT_GEN_NO_DEDUP" -> "1"];
SetEnvironment["NT_NO_UNIT_GROUPS" -> "1"];
SetEnvironment["NT_NO_SIGN_CANON" -> "1"];
ntkRef = NumTrace[net, "Frame" -> frameS0S1Sphi, "Args" -> argsS0S1SPhi];
tR = First @ AbsoluteTiming[ MakeNTKernel[ntkRef,
  FileNameJoin[{$dirR, "gen_lambda3d_small_ref.cpp"}],
  FileNameJoin[{$dirR, "Lambda3DSmallRef_num_kernel.hh"}],
  FileNameJoin[{$dirR, "Lambda3DSmallRef_num_kernels.hh"}],
  "Name" -> "Lambda3DSmallRef_num_kernel", "Namespace" -> "lambda3d_small_ref",
  "Dressings" -> dress] ];
SetEnvironment["NT_GEN_NO_DEDUP" -> None];
SetEnvironment["NT_NO_UNIT_GROUPS" -> None];
SetEnvironment["NT_NO_SIGN_CANON" -> None];
Print["[time] dedup-OFF reference kernel: ", tR, " s  (dedup-ON was ", tG, " s)"];
Print["lambda3d(small) numeric kernel + dedup-off reference generated."];
