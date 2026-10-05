#!/usr/bin/env wolframscript
(* Codegen gate for a spinor line closed on ONE head (DSL.m splitSelfTraces). The chiral condensate
   tr S(l) traces a propagator against itself: no projector closes the line, its two spinor labels
   are simply identified, e.g. ntDeltaDirac[s, s]. splitSelfTraces relabels the repeat and must insert
   the DIRAC delta; a Lorentz metric there leaves an open index and the generator aborts.

     condensate:  Mq tr[1] + Zq tr[lslash]          = 4 Mq             (the odd trace vanishes)
     product:     tr[1] tr[lslash pslash]            = 16 l.p
   compare_selftrace_num.cpp checks both values.

   Run from anywhere:  wolfram -script tests/gen/gen_selftrace_numeric.wls                          *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3;

frame = NumTracer`propFrameFT[Global`p0, Global`p, Global`l0, Global`l1, Global`cos1, Global`pm, Global`lm];
args = {Global`l0, Global`l1, Global`cos1, Global`p0, Global`p};

netA = Global`Mq[Global`l1] ntDeltaDirac[s1, s1] + Global`Zq[Global`l1] ntGamma[n1, s1, s1] ntVec[Global`lm, n1];
netB = ntDeltaDirac[s2, s2] ntGamma[n2, t1, t2] ntVec[Global`lm, n2] ntGamma[n3, t2, t1] ntVec[Global`pm, n3];

MakeNTKernel[NumTrace[netA, "Frame" -> frame, "Args" -> args, "Dressings" -> {Global`Mq, Global`Zq}],
  FileNameJoin[{$gen, "gen_selftracea_num.cpp"}],
  FileNameJoin[{$gen, "Selftracea_num_kernel.hh"}],
  FileNameJoin[{$gen, "Selftracea_num_kernels.hh"}],
  "Name" -> "Selftracea_num_kernel", "Namespace" -> "selftracea_num", "Dressings" -> {Global`Mq, Global`Zq}];

MakeNTKernel[NumTrace[netB, "Frame" -> frame, "Args" -> args, "Dressings" -> {}],
  FileNameJoin[{$gen, "gen_selftraceb_num.cpp"}],
  FileNameJoin[{$gen, "Selftraceb_num_kernel.hh"}],
  FileNameJoin[{$gen, "Selftraceb_num_kernels.hh"}],
  "Name" -> "Selftraceb_num_kernel", "Namespace" -> "selftraceb_num", "Dressings" -> {}];

Print["Dirac self-trace codegen kernels generated."];
