#!/usr/bin/env wolframscript
(* Codegen gate for FINITE-T SPATIAL VECTORS — FormTracer's `vecs[q, mu]` acronym, and above all the
   SPATIAL PSLASH `vecs[q,mu] gamma[mu,d1,d2]` that a finite-T quark line is written with.

   WHAT IS UNDER TEST. `vecs` becomes ntVec[ntSpatialVec[q], mu]: a momentum LEAF whose frame
   components are the parent's with the temporal slot zeroed (DSL.m expandSpatialVecs /
   spatialVecFrame). Nothing downstream knows about it — a spatial slash is emitted as an ordinary
   slash, and frameMask prunes the zero temporal component. That is the claim this gate has to be
   able to falsify, so the nets below are chosen to separate "spatial" from "full" everywhere:

     spat2  tr( p̸̄ (p-l)̸̄ )              = 4 [ pS.pS - pS.lS ]
            — also exercises LINEARITY: ntSpatialVec[pp - ll] must push through the difference, or
              the frame grows a third leaf whose components duplicate the other two.
     mix4   tr( p̸̄ l̸ l̸̄ p̸ )               = 4 [ (pS.lS)^2 - (pS.lS)(p.l) + (pS.pS)(lS.lS) ]
            — barred and unbarred legs INTERLEAVED. The middle term carries the FULL dot p.l while
              its partners carry spatial ones, so projecting the wrong leg (or projecting a
              neighbour along with it) moves this kernel even though the outer legs look right.
     cross  tr( p̸̄ l̸̄ ) tr( p̸ p̸ ) - 4 sps(p,l) tr( p̸ p̸ )    == 0
            — the two independent implementations of "spatial" required to agree: ntSPS, a scalar
              the frame folds into the coefficient, against the ntSpatialVec leaf that goes through
              the engine. Note it is a ZERO identity, so it is graded ABSOLUTELY against the scale
              of the other two kernels in compare_spatialvec_num.cpp — never pointwise-relative.
     leg    l̄ . P_T(p) . l̄                    = lS.lS - (lS.p)^2/(p.p)
            — the spatial vector as a PLAIN TENSOR LEG rather than a slash. This is not redundant
              with the three above: a slash is emitted as slash(Momentum{{{1.0, base}}}), which carries NO
              frameMask and reads all four components out of the env, whereas an open leg is emitted
              as lvec<id, base, MASK> and the mask comes from the FRAME while the components come
              from the frame SPEC. Only this net makes those two agree in anger — and it uses the
              loop's spatial vector, whose mask (0b0110) genuinely differs from its parent's
              (0b0111). It also puts a 1/q^2 atom on the parent momentum next to the spatial child.

   WHY A NON-ZERO CONTROL MATTERS. `vecs` support failing OPEN (a spatial vector silently treated as
   the full one) leaves every kernel finite and plausible; failing CLOSED (the leaf traced to
   nothing) makes `cross` pass for free. So spat2/mix4 are non-zero by construction and are checked
   against a closed-form oracle, and `cross` is only believed once they are large.

   FRAME. propFrameFT with a NONZERO external Matsubara component p0 and a nonzero loop l0. At
   p0 = l0 = 0 the spatial and the full vector coincide and every statement here holds for the wrong
   reason — the same trap gen_emvertex_numeric.wls had to leave the symmetric point for. *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

frame = NumTracer`propFrameFT[Global`p0, Global`p, Global`l0, Global`l1, Global`cos1,
                              Global`pp, Global`ll];
args  = {Global`p0, Global`p, Global`l0, Global`l1, Global`cos1};

(* the spatial part of a momentum, i.e. FormTracer's vecs[q, mu] leg *)
sv[q_] := ntSpatialVec[q];

(* a closed Dirac trace of the given slash momenta, tr(q̸1 q̸2 ...): one gamma per momentum, spinor
   labels chained back to the start. Each call mints its own labels so two traces can multiply. *)
tr[moms__] := Module[{qs = {moms}, n, ds, ms},
  n  = Length[qs];
  ds = Table[Unique["ds"], {n}];
  ms = Table[Unique["ms"], {n}];
  Times @@ Table[
    ntGamma[ms[[i]], ds[[i]], ds[[Mod[i, n] + 1]]] ntVec[qs[[i]], ms[[i]]], {i, n}]];

nets = <|
  (* NON-ZERO: spatial-only, and linear in a COMPOSITE spatial momentum *)
  "spat2" -> tr[sv[Global`pp], sv[Global`pp - Global`ll]],
  (* NON-ZERO: barred and unbarred legs interleaved in one four-gamma trace *)
  "mix4"  -> tr[sv[Global`pp], Global`ll, sv[Global`ll], Global`pp],
  (* ZERO: ntSPS (scalar path) against the ntSpatialVec leaf (engine path) *)
  "cross" -> tr[sv[Global`pp], sv[Global`ll]] tr[Global`pp, Global`pp] -
             4 ntSPS[Global`pp, Global`ll] tr[Global`pp, Global`pp],
  (* NON-ZERO: the spatial vector as an OPEN Lorentz leg (lvec + frameMask), not a slash *)
  "leg"   -> ntVec[sv[Global`ll], Global`a1] ntTransProj[Global`pp, Global`a1, Global`a2] *
             ntVec[sv[Global`ll], Global`a2]
|>;

emit[key_String, nm_String, ns_String] := Module[{ntk, fr, svKeys},
  ntk = NumTrace[nets[key], "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
  Print["[spatialvec] ", key, ": ", Length[ntk[[1]]["Diagrams"]], " diagram(s)"];
  (* A net that traced to NOTHING would confirm the `cross` identity with an empty kernel rather
     than by cancellation, and would make the oracle comparison vacuous for the other two. *)
  If[Length[ntk[[1]]["Diagrams"]] < 1,
    Print["FAIL: net \"", key, "\" produced NO diagrams."]; Exit[1]];
  (* THE MECHANISM ITSELF, asserted where it is cheapest to see: every spatial vector must have
     reached the frame as its own momentum with a ZERO temporal component, and must hold an env
     Base (without one, compileDirac's slash emission would not find it). *)
  fr = ntk[[1]]["Frame"];
  svKeys = Select[Keys[fr], MatchQ[#, _ntSpatialVec] &];
  If[svKeys === {},
    Print["FAIL: net \"", key, "\" has no ntSpatialVec frame entry — the spatial vector was ",
          "rewritten away before it reached the frame."]; Exit[1]];
  If[! AllTrue[svKeys, fr[#][[1]] === 0 &],
    Print["FAIL: an ntSpatialVec frame entry has a NONZERO temporal component: ",
          Select[svKeys, fr[#][[1]] =!= 0 &] -> (fr[#] & /@ svKeys)]; Exit[1]];
  (* ...and it must reach the ENV as a momentum. Not necessarily leaf-by-leaf: a composite spatial
     momentum (ntSpatialVec[pp] - ntSpatialVec[ll], what vecs[p-l,mu] becomes) is ONE env entry, so
     the test is that some env momentum mentions a spatial vector — a leaf-wise test would fail on
     exactly the linearity case spat2 exists to cover. Without an env Base compileDirac's slash
     emission would not find the momentum at all. *)
  If[! AnyTrue[Keys[ntk[[1]]["Env"]], ! FreeQ[#, _ntSpatialVec] &],
    Print["FAIL: no env momentum mentions an ntSpatialVec; env = ",
          Keys[ntk[[1]]["Env"]]]; Exit[1]];
  MakeNTKernel[ntk,
    FileNameJoin[{$gen, "gen_" <> ns <> ".cpp"}],
    FileNameJoin[{$gen, nm <> ".hh"}],
    FileNameJoin[{$gen, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> {}];
  ];

emit["spat2", "SpatialVec_spat2_kernel", "spatialvec_spat2"];
emit["mix4",  "SpatialVec_mix4_kernel",  "spatialvec_mix4"];
emit["cross", "SpatialVec_cross_kernel", "spatialvec_cross"];
emit["leg",   "SpatialVec_leg_kernel",   "spatialvec_leg"];

(* The spatial vectors must have reached BOTH emission paths — that is the whole design claim ("a
   spatial vector is an ordinary momentum, so a spatial pslash needs no new engine token"). Without
   this the gate would still pass on some future implementation that routed them elsewhere. *)
gsrc = StringJoin[(Import[#, "Text"] &) /@
   Select[FileNames["gen_spatialvec_*.cpp", $gen], FileExistsQ]];
If[! StringContainsQ[gsrc, "slash(Momentum{"],
  Print["FAIL: emitted generator source has no slash(...) — the spatial vectors did not go down ",
        "the slash path."]; Exit[1]];
If[! StringContainsQ[gsrc, "lvec<"],
  Print["FAIL: emitted generator source has no lvec<...>() — the `leg` net's spatial vector did ",
        "not go down the open-leg path, so frameMask is still untested."]; Exit[1]];
Print["[spatialvec] ok: spatial vectors emitted as slash and lvec leaves."];
Print["Finite-T spatial-vector codegen kernels generated."];
