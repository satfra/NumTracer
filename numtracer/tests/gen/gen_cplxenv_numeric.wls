#!/usr/bin/env wolframscript
(* Codegen gate for COMPLEX FRAMES at finite density: a quark leg whose Euclidean frequency is
   p0 - I muq. Two things become complex that the fill env (a `double f[]`) holds as one real slot:

     (a) the projector denominator 1/k^2 of an internal gluon line k = l - p mixing the loop with
         the quark leg: k^2 = A + i B;
     (b) a dressing atom whose value is complex, here the scalar part (q0 + I muq) of a dressed quark
         propagator numerator.

   Each complex leaf is split into a (Re, Im) pair of real slots and its powers are expanded
   binomially into the (complex) monomial coefficients, so the trace bodies stay real arithmetic.

   compare_cplxenv_num.cpp grades the kernels against a direct std::complex evaluation of the same
   integrand, at muq != 0 and at muq = 0.

   A third kernel, "Cplxcoef", pins the complex numeric COEFFICIENT of a Lorentz sum (fix: cppCx):
   the coefficient I of a sum term must be emitted as Cx{0,1}, not as CForm's `Complex(0,1)`.

   Run from anywhere:  wolfram -script tests/gen/gen_cplxenv_numeric.wls                            *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3;

(* Finite-T propagator frame with a COMPLEX external frequency. km = lm - pm is the internal gluon
   line: a frame momentum of its own so the transverse projector can carry it. *)
frame0 = NumTracer`propFrameFT[Global`p0 - I Global`muq, Global`p, Global`l0, Global`l1, Global`cos1,
                               Global`pm, Global`lm];
frame = Append[frame0, Global`km -> frame0[Global`lm] - frame0[Global`pm]];

args = {Global`l0, Global`l1, Global`cos1, Global`p0, Global`p, Global`muq};

(* (a) tr[ gamma_mu lslash gamma_nu pslash ] P^T_{mu nu}(l - p): the complex 1/k^2 *)
netA = ntGamma[m1, s1, s2] ntGamma[m3, s2, s3] ntVec[Global`lm, m3] ntGamma[m2, s3, s4] *
       ntGamma[m4, s4, s1] ntVec[Global`pm, m4] ntTransProj[Global`km, m1, m2];

(* (b) a dressed quark numerator  Zq (lslash + (l0 + I muq)) closed against (pslash + 2): the dressing
   atom (l0 + I muq) is complex and multiplies the even trace tr[1 (pslash + 2)]. The overall I puts the
   atom's imaginary half into the real kernel value. Zq is a runtime dressing, so the sum is collected, not distributed. *)
netB = I (Global`Zq[Global`l1] ntGamma[n1, t1, t2] ntVec[Global`lm, n1] +
        Global`Zq[Global`l1] (Global`l0 + I Global`muq) ntDeltaDirac[t1, t2]) *
       (ntGamma[n2, t2, t1] ntVec[Global`pm, n2] + 2 ntDeltaDirac[t2, t1]);

ntkA = NumTrace[netA, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
ntkB = NumTrace[netB, "Frame" -> frame, "Args" -> args, "Dressings" -> {Global`Zq}];

Print["[gen] complex projector denominator ..."];
MakeNTKernel[ntkA,
  FileNameJoin[{$gen, "gen_cplxenva_num.cpp"}],
  FileNameJoin[{$gen, "Cplxenva_num_kernel.hh"}],
  FileNameJoin[{$gen, "Cplxenva_num_kernels.hh"}],
  "Name" -> "Cplxenva_num_kernel", "Namespace" -> "cplxenva_num", "Dressings" -> {},
  "ComplexEndProjection" -> True, "RealOutput" -> True];

Print["[gen] complex dressing atom ..."];
MakeNTKernel[ntkB,
  FileNameJoin[{$gen, "gen_cplxenvb_num.cpp"}],
  FileNameJoin[{$gen, "Cplxenvb_num_kernel.hh"}],
  FileNameJoin[{$gen, "Cplxenvb_num_kernels.hh"}],
  "Name" -> "Cplxenvb_num_kernel", "Namespace" -> "cplxenvb_num", "Dressings" -> {Global`Zq},
  "ComplexEndProjection" -> True, "RealOutput" -> True];

(* (c) a Lorentz SUM with a complex numeric coefficient, in the vacuum frame. The overall I moves that
   coefficient's imaginary part into the real kernel value, so the grader sees it. *)
frameV = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
                             Global`q1, Global`q2, Global`q3, Global`ql];
argsV = {Global`l1, Global`cos1, Global`cos2, Global`p};
netC = I ntGamma[c1, u1, u2] ntVec[Global`ql, c1] ntGamma[c2, u2, u1] ntVec[Global`q1, c2] *
       ntVec[Global`q2, c3] (ntMetric[c3, c4] + I ntVec[Global`q3, c3] ntVec[Global`ql, c4]) ntVec[Global`q1, c4];
ntkC = NumTrace[netC, "Frame" -> frameV, "Args" -> argsV, "Dressings" -> {}];

Print["[gen] complex Lorentz-sum coefficient ..."];
MakeNTKernel[ntkC,
  FileNameJoin[{$gen, "gen_cplxcoef_num.cpp"}],
  FileNameJoin[{$gen, "Cplxcoef_num_kernel.hh"}],
  FileNameJoin[{$gen, "Cplxcoef_num_kernels.hh"}],
  "Name" -> "Cplxcoef_num_kernel", "Namespace" -> "cplxcoef_num", "Dressings" -> {},
  "ComplexEndProjection" -> True, "RealOutput" -> True];

Print["Complex-frame codegen kernels generated."];
