#!/usr/bin/env wolframscript
(* Codegen gate for TRACE FUSION of diag-dressed diagrams (CodegenKernel.m ntGroupTraces).

   A per-component dressed SU(N) line (ntSUNDiagFund) folds to a colour token multiplied into the
   integrand. Diagrams sharing their dressing coefficient AND their token are interchangeable in the
   integrand assembly, so they fuse into ONE trace group, like undressed diagrams do. A token that is
   a plain number (a pure pin, tau3/2) is folded into the dressing coefficient first.

   Four diagrams, one SU(2) flavour loop each, with different Dirac/Lorentz structures:
     runtime token  {1 -> Gu, 2 -> Gd}       -> Gu^2 + Gd^2, the same for all four: ONE trace group
     constant pin   {1 -> 1/2, 2 -> -1/2}    -> 1/2, folded into the coefficient:  ONE trace group
   Each is paired with a REFERENCE kernel spelling the flavour sum as an explicit scalar, with no
   flavour structure; compare_diagfuse_num.cpp checks the two agree pointwise.

   Run from anywhere:  wolfram -script tests/gen/gen_diagfuse_numeric.wls                           *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3;

frame = NumTracer`propFrameFT[Global`p0, Global`p, Global`l0, Global`l1, Global`cos1, Global`pm, Global`lm];
args = {Global`l0, Global`l1, Global`cos1, Global`p0, Global`p};

(* the four Dirac/Lorentz structures, each a closed spinor loop on its own labels *)
diracs[k_] := {
  ntGamma[Symbol["a" <> k], Symbol["s" <> k], Symbol["t" <> k]] ntVec[Global`lm, Symbol["a" <> k]] *
    ntGamma[Symbol["b" <> k], Symbol["t" <> k], Symbol["s" <> k]] ntVec[Global`pm, Symbol["b" <> k]],
  ntGamma[Symbol["a" <> k], Symbol["s" <> k], Symbol["t" <> k]] ntVec[Global`lm, Symbol["a" <> k]] *
    ntGamma[Symbol["b" <> k], Symbol["t" <> k], Symbol["s" <> k]] ntVec[Global`lm, Symbol["b" <> k]],
  ntGamma[Symbol["a" <> k], Symbol["s" <> k], Symbol["t" <> k]] ntVec[Global`pm, Symbol["a" <> k]] *
    ntGamma[Symbol["b" <> k], Symbol["t" <> k], Symbol["s" <> k]] ntVec[Global`pm, Symbol["b" <> k]],
  ntDeltaDirac[Symbol["s" <> k], Symbol["t" <> k]] ntDeltaDirac[Symbol["t" <> k], Symbol["s" <> k]]};

flavLoop[spec_, k_] := ntSUNDiagFund[2, Symbol["f" <> k], Symbol["g" <> k], spec] *
                       ntSUNDiagFund[2, Symbol["g" <> k], Symbol["f" <> k], spec];

(* all four diagrams carry the same dressing coefficient Z[l1] *)
build[flav_] := Total @ MapIndexed[Global`Z[Global`l1] flav[ToString[#2[[1]]]] #1 &, diracs["x"]];

runtimeSpec = {1 -> Global`Gu[Global`l1], 2 -> Global`Gd[Global`l1]};
pinSpec     = {1 -> 1/2, 2 -> -1/2};

emit[name_, net_, dress_] := MakeNTKernel[
  NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> dress],
  FileNameJoin[{$gen, "gen_" <> name <> "_num.cpp"}],
  FileNameJoin[{$gen, ToUpperCase[StringTake[name, 1]] <> StringDrop[name, 1] <> "_num_kernel.hh"}],
  FileNameJoin[{$gen, ToUpperCase[StringTake[name, 1]] <> StringDrop[name, 1] <> "_num_kernels.hh"}],
  "Name" -> ToUpperCase[StringTake[name, 1]] <> StringDrop[name, 1] <> "_num_kernel",
  "Namespace" -> name <> "_num", "Dressings" -> dress];

dressRT = {Global`Z, Global`Gu, Global`Gd};
emit["diagfusert", build[flavLoop[runtimeSpec, #]&], dressRT];
emit["diagfusertref", (Global`Gu[Global`l1]^2 + Global`Gd[Global`l1]^2) build[1&], dressRT];
emit["diagfusepin", build[flavLoop[pinSpec, #]&], {Global`Z}];
emit["diagfusepinref", (1/2) build[1&], {Global`Z}];

(* THE fusion assertion: four diagrams with one coefficient and one token are ONE trace. *)
Do[
  With[{n = StringCount[Import[FileNameJoin[{$gen, h <> "_num_kernels.hh"}], "Text"],
                        RegularExpression["inline \\S+ tr\\d+\\("]]},
    If[n =!= 1,
      Print["FAIL: ", h, " emitted ", n, " traces, expected 1: diagrams sharing their dressing ",
            "coefficient and colour token did not fuse."]; Exit[1]];
    Print["[gen] ok: ", h, " fused into 1 trace"]],
  {h, {"Diagfusert", "Diagfusepin"}}];

Print["Diag-dressed trace-fusion kernels generated."];
