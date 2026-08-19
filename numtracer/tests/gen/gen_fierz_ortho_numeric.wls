#!/usr/bin/env wolframscript
(* Projector-orthonormality gate for a FOUR-QUARK Fierz basis, through the full NumTracer pipeline.
   <P_i, V_j> = KroneckerDelta[i,j] over all n^2 pairs of FierzCompleteNf2Nc3NJL.

   This basis is the reason the fixed-Lorentz-component feature exists. It is written in a finite-T
   3+1 split, so its structures carry gamma^0 — a Lorentz index pinned to a CONCRETE component —
   which the DSL had no notion of; NumTrace aborted with `privclash` on the literal 0 and the basis
   was simply un-traceable. See expandFixedComponents in DSL.m and the gen_fixedcomp_numeric.wls
   gate, which pins the rewrite itself against analytic oracles. This gate pins the PHYSICS end:
   that a real four-quark basis traced through the rewritten path is actually orthonormal.

   Three things make it a genuinely harder test than the AqbqDirect gates:
     * FOUR external legs and TWO independent spinor lines (the earlier gates have one), so the
       multi-loop dloopsep path is exercised;
     * a TWO-TERM inner product with a MINUS sign, 2 T1[1,2,3,4] T2[2,1,4,3] - 2 T1[1,2,3,4] T2[2,3,4,1],
       which makes the contraction itself a Plus — exactly the shape of the multi-term-projector bug;
     * momentum-INDEPENDENT contact structures, so the trace must come out constant.

   As in gen_ortho_numeric.wls, the whole matrix is contracted against a distinct-prime weight matrix
   in ONE net rather than committing n^2 kernels:
       W = Sum_ij w[i,j] <P_i,V_j>   must equal   Sum_i w[i,i]
   Distinct primes mean no plausible error cancels: one wrong pair shifts W by its own weight. *)
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
SetNc[3]; SetNf[2];   (* the Nf2 basis; at the wrong Nf the projectors come out Indeterminate *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

b = "FierzCompleteNf2Nc3NJL";
n = Length[TBGetBasis[b]];
Print["[gen] basis ", b, ", n = ", n];

(* One index list per particle, INCLUDING the momentum (MakeIndexList requires it). The projector's
   quark leg closes against the vertex's ANTI-quark leg, so the WHOLE bundle moves — pinning the
   momenta while swapping only the (Dirac,colour,flavour) triples transposes the Dirac structure. *)
s1 = {Global`p1, d1, A1, F1}; s2 = {Global`p2, d2, A2, F2};
s3 = {Global`p3, d3, A3, F3}; s4 = {Global`p4, d4, A4, F4};
idxP = {s1, s2, s3, s4};
idxA = {s2, s1, s4, s3};   (* the two terms of this basis's inner product *)
idxB = {s2, s3, s4, s1};

w[i_, j_] := Prime[n (i - 1) + j];
wantDiag = Sum[w[i, i], {i, n}];
Print["[gen] expected weighted trace Sum_i w_ii = ", wantDiag];

(* Vacuous-pass guard: without the FunKit preamble TBGetProjector/TBGetVertex stay UNEVALUATED, the
   net degenerates to a symbolic expression, and NumTrace still happily returns a diagram. *)
If[! (FreeQ[TBGetProjector[b, 1, idxP], _TBGetProjector] &&
      FreeQ[TBGetVertex[b, 1, idxA], _TBGetVertex]),
  Print["FAIL: TBGetProjector/TBGetVertex did not evaluate — the basis fixture is not loaded, so ",
        "this gate would pass vacuously."]; Exit[1]];

(* REGRESSION GUARD for a fixed UPSTREAM bug (FunKit/TensorBases). TBGetProjector/TBGetVertex used
   to emit a literal, UNEVALUATED `Unique[i]` as an internal index label — `Unique` applied to the
   structure index rather than called to mint a fresh symbol. Because that label was a pure function
   of the structure index, the DIAGONAL terms (i == j) had the projector and the vertex carrying the
   byte-identical "private" dummy, and NumTrace aborted with `privclash`. Fixed upstream 2026-07-18
   (verified: no Unique survives, and the two objects' dummies are disjoint).
   It is cheap to assert and the failure mode is nasty — a consumer without a label guard silently
   contracts the two into a wrong number, and only the diagonal is affected, which is exactly where
   an orthonormality check puts its signal. *)
If[! FreeQ[{TBGetProjector[b, 1, idxP], TBGetVertex[b, 1, idxA]}, _Unique],
  Print["FAIL: TBGetProjector/TBGetVertex emitted an unevaluated Unique[...] as an index label ",
        "(upstream TensorBases regression). On the diagonal i==j the projector and vertex then ",
        "share a byte-identical 'private' dummy and the contraction is silently wrong."]; Exit[1]];

(* ONE FTerm per pair, NOT the "2 T1 T2 - 2 T1 T2" combination that TBInfo[] reports as this basis's
   inner product. That combination is already baked into TBGetVertex: printing V_j(idxA) shows it
   carries BOTH leg orderings, (2,1,4,3) and (2,3,4,1), with opposite signs and an overall 1/2 —
   i.e. the vertex IS the antisymmetrised structure. Consequently V_j(idxB) = -V_j(idxA) identically,
   and applying the inner product a second time here gave 2 F - 2 (-F) = 4 F: every diagonal entry
   came out 4x, and the gate read 9804 = 4 * 2451. That was a bug in THIS script, not in the tracer.
   (It presents as a uniform factor precisely because the distinct-prime weighting cannot hide one.) *)
net = FromFunKit @ Sum[
   w[i, j] FTerm[TBGetProjector[b, i, idxP], TBGetVertex[b, j, idxA]],
   {i, n}, {j, n}];
Print["[gen] weighted net LeafCount = ", LeafCount[net]];

(* The basis must actually USE the feature under test, or this gate silently stops covering it.
   Checked on the POST-expansion marker, not on `ntGamma[i_Integer, ...]`: `net` here is the output of
   FromFunKit, which runs expandFixedComponents, and that rewrite is exactly what turns a fixed
   component into a slash against a unit basis vector —
       gamma^i -> ntGamma[mu,d1,d2] ntVec[ntUnitVec[i],mu],  g^{i nu} -> ntVec[ntUnitVec[i],nu]
   (DSL.m). So an integer-indexed ntGamma cannot survive to this line by construction, and looking
   for one made the gate fail whatever the basis contained. ntUnitVec is the stronger assertion
   anyway: it proves the rewrite RAN, where the old spelling only proved the input had a fixed index.
   (This was masked for a long time behind an earlier failure in this script and surfaced only once
   that was fixed; the basis itself was never at fault.) *)
fixed = DeleteDuplicates @ Cases[net, NumTracer`ntUnitVec[i_] :> i, {0, Infinity}];
If[fixed === {},
  Print["FAIL: no fixed-component gamma in the net — this basis is supposed to be written in a ",
        "3+1 split; the gate is no longer covering the fixed-component path."]; Exit[1]];
Print["[gen] ok: fixed gamma components present (as unit-vector slashes): ", fixed];

frame = NumTracer`sp4Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi,
                           Global`p1, Global`p2, Global`p3, Global`p4, Global`ql];
args  = {Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi};

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[Length[ntk[[1]]["Diagrams"]] < 1,
  Print["FAIL: the weighted Fierz orthonormality net produced 0 diagrams."]; Exit[1]];
Print["[gen] traced: ", Length[ntk[[1]]["Diagrams"]], " diagram(s)"];

MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_fierz_ortho.cpp"}],
  FileNameJoin[{$gen, "FierzOrtho_kernel.hh"}],
  FileNameJoin[{$gen, "FierzOrtho_kernels.hh"}],
  "Name" -> "FierzOrtho_kernel", "Namespace" -> "fierz_ortho", "Dressings" -> {}];

Print["Four-quark Fierz orthonormality kernel generated (want ", wantDiag, ")."];
