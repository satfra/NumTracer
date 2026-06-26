#!/usr/bin/env wolframscript
(* Codegen gate for a DISCONNECTED Dirac x Lorentz diagram: a hand-built DSL net (no FunKit / no
   FORM / no dressings) where each diagram factorizes into a closed Dirac trace, an (optional) closed
   colour trace, AND an index-disjoint pure-Lorentz scalar. These three pieces MULTIPLY (they are
   independent closed scalars of ONE 1-loop diagram) — exactly the structure that tripped the
   MakeNTKernel::disconnectmix guard (the hSigL sigma-quark vertex, QCD.nb). The numeric backend now
   factors the disconnected Lorentz scalar as its own trace and multiplies it in at assembly, matching
   the dense oracle's `coeff * tr_dirac * tr_colour * tr_lorentz`. compare_discdirac_num.cpp checks
   numeric == dense over ~200k random points.

   Three diagrams over the sp3Frame kinematics (p, l1 are SCALAR magnitudes here; the 4-momenta are
   q1, q2, q3, ql — ql carries the loop angle):
     d1 = tr(q1_slash ql_slash) * (Nc colour loop) * ( q3 . P_T(ql) . q3 )   [1 Dirac + 1 Lorentz]
        = 4 (q1.ql) * Nc * ( q3^2 - (q3.ql)^2/ql^2 )
     d2 = tr(q2_slash ql_slash) * ( q1 . P_L(ql) . q1 )                       [1 Dirac + 1 Lorentz]
        = 4 (q2.ql) * (q1.ql)^2/ql^2
     d3 = tr(q1_slash ql_slash) * tr(q2_slash q3_slash)                       [2 DIRAC components]
        = 16 (q1.ql) (q2.q3)
   d3 is the hSigL structure (diagram 136): TWO disconnected Dirac traces in ONE 1-loop diagram that
   MULTIPLY. The kernel depends on the loop angle, keeps the 1/ql^2 atom, and exercises a Dirac trace
   multiplied by a transverse Lorentz scalar, by a longitudinal Lorentz scalar, and by another Dirac
   trace — with and without a constant colour factor. *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
                                   Global`q1, Global`q2, Global`q3, Global`ql];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p};

(* trailing '+' keeps the statement open across the line break (a leading '+' would end it) *)
net = (ntGamma[a1, d1, d2] ntVec[Global`q1, a1] ntGamma[a2, d2, d1] ntVec[Global`ql, a2]) *
        (ntSUNDeltaFund[3, c1, c2] ntSUNDeltaFund[3, c2, c1]) *
        (ntVec[Global`q3, i1] ntTransProj[Global`ql, i1, i2] ntVec[Global`q3, i2]) +
      (ntGamma[b1, e1, e2] ntVec[Global`q2, b1] ntGamma[b2, e2, e1] ntVec[Global`ql, b2]) *
        (ntVec[Global`q1, j1] ntLongProj[Global`ql, j1, j2] ntVec[Global`q1, j2]) +
      (ntGamma[m1, s1, s2] ntVec[Global`q1, m1] ntGamma[m2, s2, s1] ntVec[Global`ql, m2]) *
        (ntGamma[m3, t1, t2] ntVec[Global`q2, m3] ntGamma[m4, t2, t1] ntVec[Global`q3, m4]);

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];

Print["[gen] Numeric backend ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_discdirac_num.cpp"}],
  FileNameJoin[{$gen, "Discdirac_num_kernel.hh"}],
  FileNameJoin[{$gen, "Discdirac_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Discdirac_num_kernel", "Namespace" -> "discdirac_num",
  "Dressings" -> {}];

Print["[gen] Dense oracle ..."];
MakeNTKernel[ntk, FileNameJoin[{$gen, "Discdirac_num_dense_kernel.hh"}],
  "Backend" -> "Dense", "Name" -> "Discdirac_num_dense_kernel", "Dressings" -> {}];

Print["Disconnected Dirac x Lorentz codegen kernels generated."];
