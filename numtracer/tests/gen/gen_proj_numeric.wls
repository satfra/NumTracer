#!/usr/bin/env wolframscript
(* Codegen gate for the longitudinal projector: a hand-built DSL net (no FunKit / no FORM /
   no dressings) mixing ntLongProj and ntTransProj on the loop momentum, emitted via the NUMERIC
   matrix-product backend and the DENSE oracle. compare_proj_num.cpp checks numeric == dense.

   The net is two non-telescoping diagrams over the sp3Frame kinematics:
     d1 = q1 . P_L(l) . q1   (longitudinal: (q1.l)^2 / l^2, keeps the 1/l^2 atom)
     d2 = q2 . P_T(l) . q2   (transverse:  q2^2 - (q2.l)^2 / l^2)
   so the generated kernel genuinely depends on the loop angle and exercises both projector
   leaves (lproj<> / tproj<>) and the renamed projT/projL helpers in the emitted template block. *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
                                   Global`q1, Global`q2, Global`q3, Global`ql];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p};

(* trailing '+' keeps the statement open across the line break (a leading '+' would end it) *)
net = ntVec[Global`q1, i1] ntLongProj[Global`ql, i1, i2] ntVec[Global`q1, i2] +
      ntVec[Global`q2, j1] ntTransProj[Global`ql, j1, j2] ntVec[Global`q2, j2];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];

Print["[gen] Numeric backend ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_proj_num.cpp"}],
  FileNameJoin[{$gen, "Proj_num_kernel.hh"}],
  FileNameJoin[{$gen, "Proj_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Proj_num_kernel", "Namespace" -> "proj_num",
  "Dressings" -> {}];

Print["[gen] Dense oracle ..."];
MakeNTKernel[ntk, FileNameJoin[{$gen, "Proj_num_dense_kernel.hh"}],
  "Backend" -> "Dense", "Name" -> "Proj_num_dense_kernel", "Dressings" -> {}];

(* assert the new longitudinal leaf actually reached the emitted generator source (the lproj<...>()
   net-builder leaf lives in the generator .cpp / its net units; the final kernel is straight-line
   arithmetic after the generator contracts it). *)
gsrc = StringJoin[(Import[#, "Text"] &) /@
   Select[FileNames["gen_proj_num*.cpp", $gen], FileExistsQ]];
If[! StringContainsQ[gsrc, "lproj<"],
  Print["FAIL: emitted generator source does not contain an lproj<...>() leaf"]; Exit[1]];
Print["[gen] ok: generator source contains lproj<...>() and longitudinal_projector<> (dense)"];
Print["Projector codegen kernels generated."];
