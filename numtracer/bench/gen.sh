#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# bench/gen.sh — time the generation of the QCD ZA3 (three-gluon vertex) kernel
# with the FULL quark-gluon-vertex {1,4,7} tensor basis via the NUMERIC backend.
#
# Runs the generation ROUNDS (default 4) times after one warm-up round and
# reports the AVERAGE, split into the three phases of a numeric-backend build:
#   1. Mathematica frontend  — NumTracer's own work (NumTrace analyse + the
#                              C++-string codegen inside MakeNTKernel), i.e. the
#                              time spent BEFORE emitting the generator C++.
#                              The one-off FunKit derivation is done once, before
#                              the timed loop, and is NOT counted.
#   2. g++ compile           — compiling the emitted generator C++ program.
#   3. generator run         — executing that program (reduce + rebase + lower).
#
# Output kernels go to $TMPDIR so the committed tests/gen/ fixtures stay clean.
#
# Usage:  ./gen.sh [rounds]
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NT_REPO="$(cd "${SCRIPT_DIR}/../.." && pwd)"   # outer repo root (…/NumTracer)
export NT_ROUNDS="${1:-4}"

# Prefer `wolfram -script` (the form the repo's own .wls scripts use); fall back
# to a standalone `wolframscript` if that is what is installed.
if command -v wolfram >/dev/null 2>&1; then
  WL=(wolfram -script)
elif command -v wolframscript >/dev/null 2>&1; then
  WL=(wolframscript -file)
else
  echo "error: neither 'wolfram' nor 'wolframscript' found on PATH" >&2; exit 1
fi

WLS="$(mktemp --suffix=.wls)"
trap 'rm -f "$WLS"' EXIT

cat > "$WLS" <<'EOF'
(* QCD ZA3 {1,4,7} numeric-backend generation, timed in three phases. *)
$repo   = Environment["NT_REPO"];
rounds  = ToExpression[Environment["NT_ROUNDS"]];
$gen    = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
out     = $TemporaryDirectory;   (* scratch kernel outputs; keeps tests/gen clean *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];

SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];
$RecursionLimit = 1000000;

frame = NumTracer`spFrame2[Global`p, Global`l1, Global`cos1, Global`cos2,
          Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

(* --- one-off FunKit derivation (NOT timed: this is the upstream flow, not NumTracer) --- *)
fRGA3 = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3]}] // FTruncate // FPlot // FRoute // FPrint;
projA3 = FTerm[TBGetProjector["AAAClassTrans", 1, {i1, i2, i3} /. fRGA3["1-Loop"]["ExternalIndices"]]];
traceA3 = projA3 ** (fRGA3["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
spAngles = {
  Global`cosl1p1 -> Global`cos1,
  Global`cosl1p2 -> (-Global`cos1 + Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2,
  Global`cosl1p3 -> (-Global`cos1 - Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2};
flow = SPParam[traceA3 // dressingRules];
net  = FromFunKit[flow];
Print["ZA3 147 net top terms = ", If[Head[net] === Plus, Length[net], 1]];

(* Run one full generation; return {mmFrontend, compile, run} in seconds.
   The numeric backend Prints "[time] generator compile" and "[time] generator
   run" lines; we capture stdout to extract them, and back them out of the total
   MakeNTKernel time to isolate the Mathematica-side codegen work. *)
runRound[] := Module[
  {tA, ntk, capFile, strm, origOut, capText, tN, tcc, trun, mmFront},
  {tA, ntk} = AbsoluteTiming[NumTrace[net, "Frame" -> frame, "Args" -> args]];
  capFile = FileNameJoin[{$TemporaryDirectory, "ntgen_cap.txt"}];
  origOut = $Output; strm = OpenWrite[capFile]; $Output = {strm};
  tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
    FileNameJoin[{out, "ZA3_147_bench_num.cpp"}],
    FileNameJoin[{out, "ZA3_147_bench_kernel.hh"}],
    FileNameJoin[{out, "ZA3_147_bench_kernels.hh"}],
    "Backend" -> "Numeric", "Name" -> "ZA3_147_bench", "Namespace" -> "za3_147_bench",
    "Dressings" -> dress, "AngleDefs" -> spAngles] ];
  Close[strm]; $Output = origOut;
  capText = Import[capFile, "Text"];
  tcc  = ToExpression @ First @ StringCases[capText,
           "generator compile" ~~ Shortest[___] ~~ "): " ~~ x : NumberString ~~ " s" :> x, 1];
  trun = ToExpression @ First @ StringCases[capText,
           "generator run" ~~ Shortest[___] ~~ "): " ~~ x : NumberString ~~ " s" :> x, 1];
  mmFront = tA + (tN - tcc - trun);   (* NumTracer work outside g++ + generator run *)
  {mmFront, tcc, trun}];

Print["warm-up round (discarded) ..."];
runRound[];

samples = Table[
  Module[{r = runRound[]},
    Print["round ", i, ": frontend ", NumberForm[r[[1]], {6, 2}],
          " s | compile ", NumberForm[r[[2]], {6, 2}],
          " s | run ", NumberForm[r[[3]], {6, 2}], " s"];
    r],
  {i, rounds}];

mean[col_] := Mean[samples[[All, col]]];
sd[col_]   := If[Length[samples] > 1, StandardDeviation[samples[[All, col]]], 0.];

Print[""];
Print["=== QCD ZA3 {1,4,7} numeric generation — average of ", rounds, " rounds ==="];
fmt[c_, lbl_] := Print["  ", lbl, ": ",
  NumberForm[mean[c], {7, 2}], " s  (+/- ", NumberForm[sd[c], {6, 2}], ")"];
fmt[1, "1. Mathematica frontend (NumTracer)  "];
fmt[2, "2. g++ compile (generator C++)       "];
fmt[3, "3. generator run (reduce+rebase+lower)"];
Print["  ", "------------------------------------------"];
Print["  ", "total generation                     : ",
  NumberForm[mean[1] + mean[2] + mean[3], {7, 2}], " s"];
EOF

echo "Generating QCD ZA3 {1,4,7} numeric kernel, ${NT_ROUNDS} timed rounds (+1 warm-up)..."
echo "Repo: ${NT_REPO}"
"${WL[@]}" "$WLS"
