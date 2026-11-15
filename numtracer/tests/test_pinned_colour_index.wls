#!/usr/bin/env wolframscript
(* Unit gate for the pinned-component check in NumTrace: an SU(N) head whose LABEL slot holds an
   integer must abort, because labels are contracted and a component index occurring twice is summed.

   The shape that motivates it is one block of a colour-split quark propagator under a Polyakov /
   A_0 background, delta[cA, c] delta[c, cB] with a fixed integer c. Without the check the two
   deltas share a label, the label machinery reads that as a contraction, and the blocks are summed
   into one colour-blind line: the trace is wrong at A_0 != 0 and right at A_0 = 0, so every check
   without the background passes.

   The negatives matter as much as the positives. The RANK is an integer and must stay legal; a spec
   of ntSUNDiagFund / ntSUNDiagAdj carries integer component keys and is not a label slot; and the
   pinned ADJOINT path through ntSUNTPinnable (FunKitAdapter) substitutes a fresh label for the
   integer, so its output must pass.

   Run standalone:  math -noprompt -script numtracer/tests/test_pinned_colour_index.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

$nExpected = 10;
$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);

(* The check is the Cases scan in NumTrace; run it here directly so the test does not depend on a
   full trace (and so a failure points at the scan and not at something later in NumTrace).
   The level spec is {0, Infinity}: a net that IS a single head sits at level 0, and the two
   single-head cases below are exactly what a plain Infinity would miss. *)
pinnedOf[net_] := DeleteDuplicates @ Flatten @
  Cases[net, h : $ntSUNHeadPat :> Cases[Flatten[{labelsOf[h]}], _Integer], {0, Infinity}];

T[a_, i_, j_] := ntSUNT[3, a, i, j];
dF[i_, j_] := ntSUNDeltaFund[3, i, j];
dA[a_, b_] := ntSUNDeltaAdj[3, a, b];

Print["=== 1. pinned component indices are found ==="];
check["one colour block   dF[cA,1] dF[1,cB]", pinnedOf[dF[cA, 1] dF[1, cB]] === {1}];
check["all three blocks of a split propagator",
  Sort[pinnedOf[Sum[x[c] dF[cA, c] dF[c, cB], {c, 1, 3}]]] === {1, 2, 3}];
check["a pinned index on a generator  T[a,1,j]", pinnedOf[T[a, 1, j]] === {1}];
check["a pinned adjoint label         dA[1,b]", pinnedOf[dA[1, b]] === {1}];
check["hidden inside a larger net",
  pinnedOf[T[a, i, j] dF[j, k] dF[cA, 2] dF[2, cB] ntGamma[mu, d1, d2]] === {2}];

Print["=== 2. what must NOT be flagged ==="];
check["the rank is an integer and legal", pinnedOf[T[a, i, j] dF[j, i]] === {}];
check["SU(2) flavour heads, symbolic labels",
  pinnedOf[ntSUNT[2, a, i, j] ntSUNDeltaFund[2, j, i]] === {}];
(* integer keys of a dressing spec are components, but they live outside the label slots *)
check["ntSUNDiagFund spec keys are not labels",
  pinnedOf[ntSUNDiagFund[3, cA, cB, {1 -> x1, 2 -> x2, 3 -> x3}]] === {}];
check["ntSUNDiagAdj spec keys are not labels",
  pinnedOf[ntSUNDiagAdj[3, a, b, {3 -> 1}]] === {}];
(* the supported way to pin a Cartan generator: FunKitAdapter replaces the integer adjoint index by
   a fresh label and moves the component into the spec, so nothing pinned remains in a label *)
check["ntSUNTPinnable output passes",
  pinnedOf[ntSUNTPinnable[3][3, i, j]] === {}];

Print[""];
If[nRan =!= $nExpected,
  Print["FAIL: ran ", nRan, " checks, expected ", $nExpected]; Exit[1]];
If[$fail > 0, Print["FAILED: ", $fail, " check(s)"]; Exit[1]];
Print["ALL ", nRan, " CHECKS PASSED"];

End[];
