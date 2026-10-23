#!/usr/bin/env wolframscript
(* Unit gate for the colour zero-pruning in the net build: ntSingleGeneratorLoopQ, which drops a
   vanishing branch in compileColourSum and splitColourGroups and a vanishing diagram in
   ntDiagramRecords, all before the Dirac/Lorentz net build.

   tr T^a = 0. A closed fundamental index loop that runs through generators and plain Kronecker
   deltas only, and holds exactly ONE generator, makes its whole colour product vanish. Four-quark
   and Yukawa flows produce such loops in large numbers (40% of the nets of with_mesons lambda4L4),
   and before this pruning each one went through the full net build, the generator compile and the
   C++ trace only to be multiplied by a colour factor of exactly 0.

   The rule must never fire on a nonzero colour factor, so the negatives below matter as much as the
   positives: a loop with two generators, an open line, and above all a loop through a diagonal
   DRESSED delta (ntSUNDiagFund), where tr(T^a D) need not vanish.

   Run standalone:  math -noprompt -script numtracer/tests/test_colour_zero_prune.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

$nExpected = 16;
$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);

T[a_, i_, j_] := ntSUNT[3, a, i, j];
dF[i_, j_] := ntSUNDeltaFund[3, i, j];
dA[a_, b_] := ntSUNDeltaAdj[3, a, b];

Print["=== 1. vanishing: a closed fundamental loop with exactly one generator ==="];
check["tr T^a through a delta      T[a,i,j] dF[j,i]", ntSingleGeneratorLoopQ[T[a, i, j] dF[j, i]]];
check["self-traced generator       T[a,i,i]", ntSingleGeneratorLoopQ[T[a, i, i]]];
check["loop closed by a delta chain T[a,i,j] dF[j,k] dF[k,i]", ntSingleGeneratorLoopQ[T[a, i, j] dF[j, k] dF[k, i]]];
check["tr T^a tr T^b delta^ab (the lambda4L shape)",
  ntSingleGeneratorLoopQ[T[a, i, j] dF[j, i] T[b, k, l] dF[l, k] dA[a, b]]];
check["one vanishing loop beside a nonvanishing one",
  ntSingleGeneratorLoopQ[T[a, i, j] dF[j, i] T[b, k, l] T[b, l, k]]];
check["flavour SU(2)               ntSUNT[2,..] ntSUNDeltaFund[2,..]",
  ntSingleGeneratorLoopQ[ntSUNT[2, a, i, j] ntSUNDeltaFund[2, j, i]]];
check["a number and a delta^2 beside it do not hide it",
  ntSingleGeneratorLoopQ[3 T[a, i, j] dF[j, i] dF[k, l]^2]];

Print[""];
Print["=== 2. must NOT fire ==="];
check["tr(T^a T^b) delta^ab = C_F N is not zero", ! ntSingleGeneratorLoopQ[T[a, i, j] T[b, j, i] dA[a, b]]];
check["an open line                T[a,i,j] dF[j,k]", ! ntSingleGeneratorLoopQ[T[a, i, j] dF[j, k]]];
check["a dressed diagonal delta in the loop: tr(T^a D) need not vanish",
  ! ntSingleGeneratorLoopQ[T[a, i, j] ntSUNDiagFund[3, j, i, {1 -> 1}]]];
check["no generator at all         dF[i,j] dF[j,i]", ! ntSingleGeneratorLoopQ[dF[i, j] dF[j, i]]];
check["a pure delta loop as a power dF[i,j]^2", ! ntSingleGeneratorLoopQ[dF[i, j]^2]];
check["f^abc T^a T^b T^c in one loop (three generators)",
  ! ntSingleGeneratorLoopQ[ntSUNf[3, a, b, c] T[a, i, j] T[b, j, k] T[c, k, i]]];

Print[""];
Print["=== 3. compileColourSum drops exactly the vanishing branches ==="];
ids = <|a -> 0, b -> 1, i -> 2, j -> 3, k -> 4, l -> 5|>;
With[{br = compileColourSum[T[a, i, j] dF[j, i] T[a, k, l] dF[l, k] + 2 T[a, i, j] T[a, j, i], ids]},
  check["a sum of a vanishing and a nonvanishing branch keeps one branch", Length[br] === 1];
  check["... and it is the nonvanishing one, with its scalar", Length[br] === 1 && br[[1, 2]] === 2]];
check["an entirely vanishing sum leaves no branch",
  compileColourSum[T[a, i, j] dF[j, i] T[a, k, l] dF[l, k], ids] === {}];

Print[""];
If[nRan =!= $nExpected,
  Print["FAIL: ran ", nRan, " checks, expected ", $nExpected]; Exit[1]];
If[$fail > 0, Print["FAILED: ", $fail, " check(s)"]; Exit[1]];
Print["ALL ", nRan, " CHECKS PASSED"];

End[];
