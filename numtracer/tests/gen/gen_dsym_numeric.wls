#!/usr/bin/env wolframscript
(* Codegen gate for the SU(N) totally symmetric tensor ntSUNd (d^{abc} = 2 tr({T^a,T^b} T^c), the
   SUNFacKind::D factor of network/sun_net.hpp). FunKit has no d token, so the head is written by hand.
   With tr(T^a T^b) = delta^{ab}/2 in SU(3):

     dd:    Zq tr[lslash pslash] d^{abc} d^{abc}               = 4 l.p Zq (N^2-4)(N^2-1)/N = 160/3 l.p Zq
     dttt:  Mq tr[1] d^{abc} tr(T^a T^b T^c)                   = 4 Mq d.d/4                = 40/3 Mq
     sum:   tr[lslash pslash] d^{abc} (d^{abc} + f^{abc})      = 160/3 l.p                 (f.d = 0)
     diag:  tr[1] d^{abc} d^{a'b'c} diag_adj(a,a') diag_adj(b,b')
            with only lambda8 kept on a and only lambda4 on b    = 4 ZA8 ZA4 (d^{448})^2     = ZA8 ZA4 / 3

   dd and dttt check the emission and the d table (dttt also against the generator traces), sum the
   colour-Plus lowering, and diag the per-component dressed fold. diag is component-sensitive: a
   permuted d table or component order would not give d^{448}^2 = 1/12. A diag spec index is 1-based in
   NumTracer's generator order (off-diagonal first, Cartan last), so 3 is lambda4 and 8 is lambda8.
   compare_dsym_num.cpp checks all four.

   Run from anywhere:  wolfram -script tests/gen/gen_dsym_numeric.wls                                *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3;

frame = NumTracer`propFrameFT[Global`p0, Global`p, Global`l0, Global`l1, Global`cos1, Global`pm, Global`lm];
args = {Global`l0, Global`l1, Global`cos1, Global`p0, Global`p};

trLP[n1_, n2_, s1_, s2_] := ntGamma[n1, s1, s2] ntVec[Global`lm, n1] ntGamma[n2, s2, s1] ntVec[Global`pm, n2];

netDD   = Global`Zq[Global`l1] trLP[n1, n2, s1, s2] ntSUNd[3, a, b, c] ntSUNd[3, a, b, c];
(* parenthesised: a script ends a statement at a line break when the line is already complete *)
netDTTT = (Global`Mq[Global`l1] ntDeltaDirac[s1, s2] ntDeltaDirac[s2, s1] ntSUNd[3, a, b, c]
            ntSUNT[3, a, i1, i2] ntSUNT[3, b, i2, i3] ntSUNT[3, c, i3, i1]);
netSum  = trLP[n1, n2, s1, s2] ntSUNd[3, a, b, c] (ntSUNd[3, a, b, c] + ntSUNf[3, a, b, c]);
netDiag = (ntDeltaDirac[s1, s2] ntDeltaDirac[s2, s1] ntSUNd[3, a, b, c] ntSUNd[3, a2, b2, c]
            ntSUNDiagAdj[3, a, a2, {8 -> Global`ZA8[Global`l1]}] ntSUNDiagAdj[3, b, b2, {3 -> Global`ZA4[Global`l1]}]);

gen[net_, tag_String, dressings_List] :=
  MakeNTKernel[NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> dressings],
    FileNameJoin[{$gen, "gen_dsym" <> tag <> "_num.cpp"}],
    FileNameJoin[{$gen, "Dsym" <> tag <> "_num_kernel.hh"}],
    FileNameJoin[{$gen, "Dsym" <> tag <> "_num_kernels.hh"}],
    "Name" -> "Dsym" <> tag <> "_num_kernel", "Namespace" -> "dsym" <> tag <> "_num", "Dressings" -> dressings];

gen[netDD, "dd", {Global`Zq}];
gen[netDTTT, "dttt", {Global`Mq}];
gen[netSum, "sum", {}];
gen[netDiag, "diag", {Global`ZA8, Global`ZA4}];

Print["symmetric d^{abc} codegen kernels generated."];
