#!/usr/bin/env wolframscript
(* QCD generator: emit the gluon self-energy (ZA) flow. Unlike the Yang-Mills YM_ZA, the QCD
   field space (gluon+ghost+quark) adds a QUARK LOOP, so this flow mixes a Dirac-trace quark
   loop with the gluon/ghost Lorentz+colour diagrams. Nf=2 is baked in at generation time
   (qcd_setup.m calls SetNf[2]); the kernel carries no Nf parameter. 1 loop angle.
   Produces the invariant-basis kernel (production) and a dense baseline.
   Run from numtracer/:   wolfram -script tests/gen/gen_qcd_za.wls                            *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];  (* tests: emit against the shim + DiFfRG namespace *)

(* symmetric-point frame, 1 angle: p1 along axis 0, loop l1 at angle cos1 *)
frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
args  = {Global`l1, Global`cos1, Global`p, Global`k};
(* full QCD dressing set (matches the FORM ZA oracle parameter list) *)
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

fRGAA = FTakeDerivatives[WetterichEquation, {A[i1], A[i2]}] // FTruncate // FPlot // FRoute // FPrint;
(* the FORM oracle generates FlowAA/p^2; the 1/p^2 enters here as the FTerm scalar coeff *)
trace = FTerm[1/sp[p1, p1], TBGetProjector["AA", 1, {i1, i2} /. fRGAA["1-Loop"]["ExternalIndices"]]] **
        (fRGAA["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]) //. p2 -> -p1;
flow = PropParam[trace // dressingRules];     (* vacuum 1-angle scalar reduction -> {l1,p,cos1} *)
net  = FromFunKit[flow];                        (* contractFlavour folds the quark loop -> Nf=2 *)
Print["ZA net: top terms = ", If[Head[net] === Plus, Length[net], 1]];
Print["Nf baked as integer = ", Global`Nf, "  (SetNf[2] => folds into the quark-loop coefficient)"];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
Print["[time] NumTrace (analyse): ", tA, " s"];

(* The build-time generation path is the NUMERIC backend — see gen_qcd_za_numeric.wls
   (ZA_num_kernel). This generator emits only the ET + dense cross-check kernels. *)

(* dense baseline (ZA is the small 1-angle flow; dense always compiles) *)
tD = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "ZA_dense_kernel.hh"}],
  "Name" -> "ZA_dense_kernel", "Dressings" -> dress, "Backend" -> "Dense"] ];
Print["[time] MakeNTKernel (Dense): ", tD, " s"];

(* symbolic ET kernel (default Backend). ZA is a 1-angle propagator (like Zq, whose ET kernel
   compiles); the QCD quark+gluon+ghost loops may still blow the symbolic type up at C++ compile
   time — emit it and let the build decide whether it is feasible. *)
tE = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "ZA_gen_kernel.hh"}],
  "Name" -> "ZA_gen_kernel", "Dressings" -> dress] ];
Print["[time] MakeNTKernel (ET): ", tE, " s"];
Print["ZA kernels generated."];
