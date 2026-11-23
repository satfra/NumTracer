#!/usr/bin/env wolframscript
(* ntPhiEvenQ (CodegenKernel.m): the proof behind DiFfRG's `phi_even` trait. phi -> -phi must be a reflection of the
   frame that fixes every phi-free vector, the diagrams must carry no Levi-Civita tensor or gamma5, and the scalar body
   must hold phi only in Cos and even powers of Sin. Runs in a plain ctest (no generation, no C++). *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

pe = NumTracer`Private`ntPhiEvenQ;
fail = False;
check[label_, got_, want_] := If[got === want,
  Print["  ok   ", label],
  (Print["FAIL: ", label, "\n  got      ", InputForm[got], "\n  expected ", InputForm[want]]; fail = True)];

(* the W5 frame: soft legs in (e0, e1), the loop leg l in (e0, e1, e2), the loop momentum with three angles *)
loop = l1 {cos1, Sqrt[1 - cos1^2] cos2, Sqrt[(1 - cos1^2) (1 - cos2^2)] Cos[phi], Sqrt[(1 - cos1^2) (1 - cos2^2)] Sin[phi]};
frame = <|p1 -> S0 {1, 0, 0, 0}, p2 -> S0 {-1/2, Sqrt[3]/2, 0, 0}, p3 -> S0 {-1/2, -Sqrt[3]/2, 0, 0},
  p4 -> lm {Sqrt[1 - uu] Cos[psi], Sqrt[1 - uu] Sin[psi], Sqrt[uu], 0}, l1 -> loop|>;
diags = {ntVec[l1, mu] ntVec[p4, mu] ntTransProj[l1 - p1, nu, rho]};
body = RB[k^2, l1^2] / (l1^2 + 2 l1 lm Sqrt[1 - cos1^2] cos2 Cos[phi] + lm^2)^2 ZA[Sqrt[l1^2 Sin[phi]^2 + 1]];

check["the W frame, covariant diagrams, an even body", pe["phi", frame, <||>, diags, body], True];
check["no PhiEvenVar", pe[None, frame, <||>, diags, body], False];
check["a name that is no symbol of the flow", pe["chi", frame, <||>, diags, body], False];
check["l tilted into the reflected direction",
  pe["phi", Append[frame, p4 -> lm {Sqrt[1 - uu] Cos[psi], Sqrt[1 - uu] Sin[psi], Sqrt[uu/2], Sqrt[uu/2]}], <||>, diags, body], False];
check["gamma5 in a diagram", pe["phi", frame, <||>, {ntGamma5[d1, d2] diags[[1]]}, body], False];
check["a Levi-Civita tensor in a diagram", pe["phi", frame, <||>, {ntEpsilon[a, b, c, d] diags[[1]]}, body], False];
check["a Levi-Civita tensor written in NumTracer's context",
  pe["phi", frame, <||>, {NumTracer`Private`ntEpsilon[a, b, c, d] diags[[1]]}, body], False];
check["a fixed Cartesian axis", pe["phi", frame, <||>, {ntVec[ntUnitVec[3], mu] diags[[1]]}, body], False];
check["a hand-written pinned index", pe["phi", frame, <||>, {ntMetric[3, mu] diags[[1]]}, body], False];
check["Sin[phi] at an odd power in the body", pe["phi", frame, <||>, diags, body Sin[phi]], False];
check["a frame component of mixed parity",
  pe["phi", Append[frame, l1 -> ReplacePart[loop, 3 -> loop[[3]] + l1 Sin[phi]]], <||>, diags, body], False];
(* a symbol def that hides an odd Sin[phi] in the body *)
check["an odd Sin[phi] behind a symbol def", pe["phi", frame, <|ntS -> Sin[phi]|>, diags, body ntS], False];
check["an even one behind a symbol def", pe["phi", frame, <|ntS -> Sin[phi]^2|>, diags, body ntS], True];

If[fail, Print["test_phi_even: FAILED"]; Exit[1], Print["test_phi_even: all passed"]; Exit[0]];
