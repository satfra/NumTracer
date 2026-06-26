#!/usr/bin/env wolframscript
(* Codegen gate for FINITE-TEMPERATURE support: a hand-built DSL net (no FunKit / no FORM /
   no dressings) over the finite-T frame propFrameFT, exercising the new front-end features
   end-to-end through BOTH the numeric matrix-product backend and the dense oracle:

     - ntElectricProj / ntMagneticProj  (the finite-T electric/magnetic projectors; eproj/mproj
       leaves + electric_projector/magnetic_projector dense, with the spatial 1/|q⃗|² env slot),
     - ntTransProj at finite T  (still valid with a distinct temporal component),
     - ntSPS  (the spatial scalar product, resolved by the frame as the components-1..3 dot),
     - ntVec[q, 0]  (a scalar momentum component — the temporal/Matsubara component).

   The frame is the finite-T quark-self-energy kinematics: the external pp carries an independent
   temporal component p0 (e.g. the lowest fermionic Matsubara frequency πT) and a spatial magnitude
   p; the loop ll carries an independent temporal l0 and a spatial 3-vector (magnitude l1, angle
   cos1). Because ll's temporal component is not ∝ l1, the numeric backend takes the general
   polyFrameSpec path automatically (unitLoopOkQ is False). compare_ftproj_num.cpp checks
   numeric == dense over random kinematics. *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

frame = NumTracer`propFrameFT[Global`p0, Global`p, Global`l0, Global`l1, Global`cos1,
                               Global`pp, Global`ll];
args  = {Global`p0, Global`p, Global`l0, Global`l1, Global`cos1};

(* three non-telescoping diagrams (trailing '+' keeps the statement open across the break):
     d1 = pp . P_E(ll) . pp                                    (electric:  1/ll² and 1/|ll⃗|² atoms)
     d2 = pp . P_M(ll) . pp                                    (magnetic:  1/|ll⃗|² atom, temporal rows 0)
     d3 = (ll⃗·pp⃗)(ll_0) · pp . P_T(ll) . pp                   (scalar sps & vec[ll,0] × transverse) *)
net = ntVec[Global`pp, i1] ntElectricProj[Global`ll, i1, i2] ntVec[Global`pp, i2] +
      ntVec[Global`pp, j1] ntMagneticProj[Global`ll, j1, j2] ntVec[Global`pp, j2] +
      ntSPS[Global`pp, Global`ll] ntVec[Global`ll, 0] *
        ntVec[Global`pp, k1] ntTransProj[Global`ll, k1, k2] ntVec[Global`pp, k2];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];

Print["[gen] Numeric backend ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_ftproj_num.cpp"}],
  FileNameJoin[{$gen, "FTProj_num_kernel.hh"}],
  FileNameJoin[{$gen, "FTProj_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "FTProj_num_kernel", "Namespace" -> "ftproj_num",
  "Dressings" -> {}];

Print["[gen] Dense oracle ..."];
MakeNTKernel[ntk, FileNameJoin[{$gen, "FTProj_num_dense_kernel.hh"}],
  "Backend" -> "Dense", "Name" -> "FTProj_num_dense_kernel", "Dressings" -> {}];

(* assert the new finite-T projector leaves reached the emitted generator source (the eproj/mproj
   net-builder leaves live in the generator .cpp / its net units; the final kernel is straight-line
   arithmetic after the generator contracts it). *)
gsrc = StringJoin[(Import[#, "Text"] &) /@
   Select[FileNames["gen_ftproj_num*.cpp", $gen], FileExistsQ]];
If[! (StringContainsQ[gsrc, "eproj<"] && StringContainsQ[gsrc, "mproj<"]),
  Print["FAIL: emitted generator source missing eproj<...>() / mproj<...>() leaves"]; Exit[1]];
dsrc = Import[FileNameJoin[{$gen, "FTProj_num_dense_kernel.hh"}], "Text"];
If[! (StringContainsQ[dsrc, "electric_projector<"] && StringContainsQ[dsrc, "magnetic_projector<"]),
  Print["FAIL: dense kernel missing electric_projector<> / magnetic_projector<>"]; Exit[1]];
Print["[gen] ok: eproj/mproj (numeric) and electric/magnetic_projector (dense) present."];
Print["Finite-T projector codegen kernels generated."];
