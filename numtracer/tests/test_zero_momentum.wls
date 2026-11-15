#!/usr/bin/env wolframscript
(* Unit gate for resolveComponents (Frames.m) on a momentum that cancels identically.

   A leg sum can vanish in the chosen kinematics -- the t = (p1 + p4)^2 channel of a four-point
   function at a density-channel configuration, or any symmetric point with cancelling legs. Simplify
   then gives the NUMBER 0 rather than a list of components, and generation used to abort. The zero
   four-vector is what that means, and every caller of resolveComponents handles it (ntSP, ntSPS,
   ntVec, the frame environment), so it resolves to {0,0,0,0}.

   The negative matters: a momentum the frame does not describe is still an error, and the message
   must name what it resolved to.

   Run standalone:  math -noprompt -script numtracer/tests/test_zero_momentum.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

$nExpected = 7;
$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);

(* p1 + p2 + p3 = 0 at the symmetric point, so p1 + p2 = -p3 and p1 + p2 + p3 cancels. *)
frame = <|p1 -> S0 {1, 0, 0, 0}, p2 -> S0 {-1/2, Sqrt[3]/2, 0, 0}, p3 -> S0 {-1/2, -Sqrt[3]/2, 0, 0},
          l1 -> l1 {c1, Sqrt[1 - c1^2], 0, 0}|>;

Print["=== 1. a cancelling momentum is the zero four-vector ==="];
check["the literal number 0", resolveComponents[0, frame] === {0, 0, 0, 0}];
check["a leg sum that cancels   p1 + p2 + p3", resolveComponents[p1 + p2 + p3, frame] === {0, 0, 0, 0}];
check["a difference that cancels  l1 - l1", resolveComponents[l1 - l1, frame] === {0, 0, 0, 0}];
check["machine zero             0.", resolveComponents[0., frame] === {0, 0, 0, 0}];

Print["=== 2. the callers give zero on it, as they must ==="];
check["scalar product with a zero vector",
  resolveComponents[p1 + p2 + p3, frame] . resolveComponents[p1, frame] === 0];
check["a component of the zero vector", resolveComponents[p1 + p2 + p3, frame][[2]] === 0];

Print["=== 3. a momentum the frame does not describe still aborts ==="];
(* CheckAbort, not Check: the failure is an Abort[] after the message, and an uncaught one would
   take this script with it. *)
check["an unknown symbol aborts",
  TrueQ @ Quiet[CheckAbort[resolveComponents[qUnknown, frame]; False, True],
    resolveComponents::zeroMomentum]];

Print[""];
If[nRan =!= $nExpected,
  Print["FAIL: ran ", nRan, " checks, expected ", $nExpected]; Exit[1]];
If[$fail > 0, Print["FAILED: ", $fail, " check(s)"]; Exit[1]];
Print["ALL ", nRan, " CHECKS PASSED"];

End[];
