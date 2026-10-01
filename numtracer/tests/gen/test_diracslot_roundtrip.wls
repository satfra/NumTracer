#!/usr/bin/env wolframscript
(* Stage 4 R2: round-trip test for the general Dirac-slot decompose/expand (any open-leg count).
   Synthetic k=1 and k=2 vertex sums (no FunKit needed) — decompose then expand must equal Expand[sum]. *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

(* pull the private symbols under test *)
dec  = NumTracer`Private`diracSlotDecompose;
(* expand a slot back into its sum: each option is {scalar, structure} *)
exp[NumTracer`Private`ntDiracSlot[opts_, __]] := Total[Times @@@ opts];
sumQ = NumTracer`Private`diracSlotSumQ;
openL = NumTracer`Private`openLorentzOf;
nDS  = NumTracer`Private`ntDiracSlot;   (* ntDiracSlot lives in the Private context, like ntDressedNum *)

ng   = NumTracer`ntGamma; nd = NumTracer`ntDeltaDirac; nv = NumTracer`ntVec;
nmet = NumTracer`ntMetric; nT = NumTracer`ntSUNT; nProj = NumTracer`ntTransProj;

osp = NumTracer`Private`openSpinorOf;
fails = 0;
check[name_, sum_] := Module[{d, r, lhs, rhs, terms},
  Print["\n== ", name, " =="];
  terms = List @@ Expand[sum];
  Print["  per-term openLorentz = ", openL /@ terms];
  Print["  per-term openSpinor  = ", osp /@ terms];
  Print["  diracSlotSumQ = ", TrueQ[sumQ[sum]]];
  d = dec[sum];
  If[d === $Failed, Print["  decompose = $Failed  FAIL"]; fails++; Return[]];
  Print["  slot legs = ", Cases[d, nDS[_, _, _, legs_] :> legs, {0, Infinity}]];
  Print["  #options  = ", Cases[d, nDS[o_, _, _, _] :> Length[o], {0, Infinity}]];
  r = d /. s_nDS :> exp[s];
  lhs = Expand[sum]; rhs = Expand[r];
  If[Simplify[lhs - rhs] === 0,
     Print["  round-trip: expand@decompose == Expand[sum]  ok"],
     Print["  round-trip MISMATCH:\n   lhs-rhs = ", Simplify[lhs - rhs], "  FAIL"]; fails++]];

(* common factors: a colour generator and a scalar denominator, shared by every structure. *)
common = nT[3, a, i, j] * denom;

(* ---- k=1 (single-gluon vertex, AqbqDirect147-like): open leg = mu, closed by a projector ----
   T1 = γ^μ ; T4 = (p̸1)·γ^ρ with the projector on ρ (open) ; each dressing differs. *)
v1 = Z1 * common * ng[mu, d1, d3] * nProj[p1, ext, mu] +
     Z4 * common * ng[a1, d1, d2] * nv[p1, a1] * ng[mu, d2, d3] * nProj[p1, ext, mu] +
     Z7 * common * ng[a2, d1, d2] * nv[l1, a2] * ng[mu, d2, d3] * nProj[p1, ext, mu];
(* open leg = ext (the projector's external axis); μ,ρ are contracted with the projector/slashes. *)
check["k=1 single-gluon vertex (open leg via projector)", v1];

(* ---- k=2 (two-gluon vertex, AAqbq-like): open legs = {mu, nu} ----
   γ^μγ^ν ; g^{μν}·δ (both legs on a metric) ; γ^μ·(γ·q) with a p^ν vector (mixed). *)
v2 = W1 * common * ng[mu, d1, d2] * ng[nu, d2, d3] +
     W2 * common * nd[d1, d3] * nmet[mu, nu] +
     W3 * common * ng[mu, d1, d2] * ng[rho, d2, d3] * nv[q, rho] * nv[qq, nu];
check["k=2 two-gluon vertex (Dirac / metric / mixed legs)", v2];

(* ---- R3: the collection gate keeps the vertex eager and rewrites it to ntDiracSlot ---- *)
collQ = NumTracer`Private`collectibleDiracSumQ;
rew   = NumTracer`Private`rewriteDressedNums;
(* a DRESSED vertex sum is collected only under vertex collection (NT_VERTEX_COLLECT), so turn it on *)
gate[name_, sum_] := Module[{c, r},
  Block[{NumTracer`Private`$ntDressCollect = True, NumTracer`Private`$ntVertexCollect = True},
    c = TrueQ[collQ[sum]];
    r = rew[{sum}]];
  Module[{hasSlot = ! FreeQ[r, nDS]},
    Print["  ", name, ": collectibleDiracSumQ=", c, "  rewrite→ntDiracSlot=", hasSlot,
          "  ", If[c && hasSlot, "ok", "FAIL"]];
    If[! (c && hasSlot), fails++]]];
Print["\n== R3: collection gate (collectible + rewrite to ntDiracSlot) =="];
gate["k=1 vertex", v1];
gate["k=2 vertex", v2];

Print["\n", If[fails == 0, "ALL ROUND-TRIP TESTS PASSED", ToString[fails] <> " FAILED"]];
Exit[If[fails == 0, 0, 1]];
