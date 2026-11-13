#!/usr/bin/env wolframscript
(* canonicalizeEqualMomenta (DSL.m): momenta that are equal as frame vectors are rewritten to one
   spelling, so that they share one env Base and Inv slot. Here p1 + p2 + p3 = 0, so l1 - p1 and
   l1 + p2 + p3 are the same vector. Runs in a plain ctest (no generation, no C++). *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3;

cmf = NumTracer`Private`canonicalizeEqualMomenta;
fail = False;
check[label_, got_, want_] := If[got === want,
  Print["  ok   ", label],
  (Print["FAIL: ", label, "\n  got      ", InputForm[got], "\n  expected ", InputForm[want]]; fail = True)];

frame = <|p1 -> S0 {1, 0, 0, 0}, p2 -> S0 {-1/2, Sqrt[3]/2, 0, 0}, p3 -> S0 {-1/2, -Sqrt[3]/2, 0, 0},
          l1 -> l1 {c1, s1, 0, 0}|>;

(* --- A. The two spellings of one vector become one; the first met (in Cases order) is kept. ------- *)
check["l1 + p2 + p3 and l1 - p1 share one spelling",
  cmf[ntTransProj[l1 + p2 + p3, mu, nu] ntVec[l1 - p1, mu] ntVec[p2, nu], frame],
  ntTransProj[l1 + p2 + p3, mu, nu] ntVec[l1 + p2 + p3, mu] ntVec[p2, nu]];

(* --- B. The negative of a vector: a linear head takes its sign, a projector does not. Which of the
   two spellings is kept depends on the factor order, so the check holds for either. --------------- *)
With[{out = cmf[ntVec[p1 - l1, mu] ntVec[l1 + p2 + p3, nu], frame]},
  check["p1 - l1 = -(l1 + p2 + p3): a vec flips sign",
    MemberQ[{-ntVec[#, mu] ntVec[#, nu]} & /@ {p1 - l1, l1 + p2 + p3}, {out}], True]];
With[{out = cmf[ntTransProj[p1 - l1, mu, nu] ntTransProj[l1 + p2 + p3, nu, rho], frame]},
  check["p1 - l1 = -(l1 + p2 + p3): a projector does not",
    MemberQ[{ntTransProj[#, mu, nu] ntTransProj[#, nu, rho]} & /@ {p1 - l1, l1 + p2 + p3}, {out}], True]];

(* --- C. Nothing to merge, and a momentum the frame does not resolve: unchanged. ------------------ *)
With[{net = ntTransProj[l1 - p1, mu, nu] ntVec[p2, mu] ntVec[q, nu]},
  check["distinct vectors and an unresolved momentum: unchanged", cmf[net, frame], net]];

(* --- D. Sigma slash legs carry {coeff, q} pairs; both legs are rewritten. ------------------------ *)
check["both slash legs of a sigma",
  cmf[ntVec[l1 - p1, mu] ntSigma[{"slash", {{1, l1 + p2 + p3}}}, {"slash", {{2, p1 - l1}}}, d1, d2], frame],
  ntVec[l1 - p1, mu] ntSigma[{"slash", {{1, l1 - p1}}}, {"slash", {{-2, l1 - p1}}}, d1, d2]];

If[fail, Print["equal momenta: FAILED"]; Exit[1]];
Print["equal momenta: all checks passed"];
