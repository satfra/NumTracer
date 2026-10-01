#!/usr/bin/env wolframscript
(* Spinor-loop ORIENTATION guard (Codegen.m: orderDiracFacs, DSL.m: diracSlotDecompose).

   A closed spinor loop is emitted as a left-to-right 4x4 MATRIX PRODUCT tr(M1..Mn) — that collapse
   is what makes phase A affordable — and a matrix product has an ORDER, which orderDiracFacs
   recovers by walking the loop. The walk is UNDIRECTED, so it always closes: it just follows edges.
   That means the pre-existing "a spinor loop did not close" guard cannot see a loop whose fermion
   arrow is inconsistent (two `in` slots meeting, as at an anomalous qq / q̄q̄ diquark vertex). Such a
   loop needs one SEGMENT transposed. That is now REPRESENTABLE: a factor the walk entered against
   its declared (diracIn,diracOut) order is wrapped in `ntTransposed`, and the engine multiplies its
   transpose. So this test pins the MARKING, not a refusal — the refusal that remains is only for a
   collected slot, whose transpose is not implemented.

   Asserting on Length[] alone would be a FALSE GREEN: the wrapper preserves the token count, so a
   walk that marked nothing would still pass. Every case below checks the flags themselves.

   Runs in a PLAIN ctest (no generation, no C++), because codegen_regen is opt-in and in practice is
   not run while editing. *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3;

odf  = NumTracer`Private`orderDiracFacs;
fail = False;
check[label_, got_, want_] := If[got === want,
  Print["  ok   ", label, " -> ", got],
  (Print["FAIL: ", label, " -> ", got, ", expected ", want]; fail = True)];
(* the transpose marking of each returned token, in order *)
flags[out_] := MatchQ[#, _NumTracer`Private`ntTransposed] & /@ out;

(* --- A. A consistently oriented loop: every token forward, nothing marked. ------------------- *)
check["forward loop: nothing marked",
  CheckAbort[flags[odf[{ntGamma[a1, d1, d2], ntGamma[a2, d2, d3], ntGamma[a3, d3, d1]}]], "ABORTED"],
  {False, False, False}];

(* --- B. A UNIFORMLY reversed loop is legal and must NOT abort. ---------------------------------
   This is not hypothetical: the walk seeds by entering facs[[1]] on its `diracIn`, and when that
   factor is the symmetric external projector ntDeltaDirac — whose two labels carry no orientation —
   the arbitrary label order decides which way round the cycle goes. The sigL / etaPiL meson
   closures traverse EVERY factor backwards for exactly this reason. It costs nothing, because
   reversing a closed gamma loop is the charge-conjugation identity: tr(Mn..M1) = tr(M1^T..Mn^T) and
   gamma_mu^T = C^-1 (-gamma_mu) C, so a reversal is worth (-1)^(#gamma) — and any trace with an odd
   gamma count already vanishes by the parity rule applied elsewhere in Codegen.m. (gamma5 and the
   identity transpose with a PLUS sign and do not enter that count.) So the engine has ALWAYS relied
   on this identity implicitly; a guard keyed on direction rather than on mixture would break every
   meson flow in the tree. Two tokens: the delta connector carries none. *)
check["uniformly reversed loop: every token marked",
  CheckAbort[flags[odf[{ntDeltaDirac[d1, d2], ntGamma[a2, d3, d2], ntGamma[a3, d1, d3]}]], "ABORTED"],
  {True, True}];

(* --- C. A MIXED loop — the arrow structure of a diagram with an anomalous qq vertex and its
   q̄q̄ conjugate: two IN ends meet at d2, two OUT ends at d4. It is no longer refused; the two
   backward factors are marked and the engine transposes them. This is the diquark case. -------- *)
check["anomalous (mixed) loop marks exactly the backward factors",
  CheckAbort[flags[odf[{ntGamma[a1, d1, d2], ntGamma[a2, d3, d2],
                        ntGamma[a3, d3, d4], ntGamma[a4, d1, d4]}]], "ABORTED"],
  {False, True, False, True}];

(* --- D. A COLLECTED token traversed backwards is marked too. It emits `dtrslot(k)` rather than
   `dtr(...)`: a slot is a sum of chains, and (A1..An)^T = An^T..A1^T, so the reversal happens at the
   C++ dress_enumerate splice. The transpose rides the slot REFERENCE, not the slot, so the
   $dslCache key is unchanged. ------------------------------------------------------------------ *)
(* ntDressedNum has no ::usage, so it lives in NumTracer`Private` and must be qualified — an
   unqualified spelling silently makes a DIFFERENT symbol that matches no head pattern, and the walk
   then reports a loop that did not close. *)
check["reversed collected slot is marked",
  CheckAbort[flags[odf[{ntDeltaDirac[e1, e2], NumTracer`Private`ntDressedNum[{{1, {"ident"}}}, e3, e2],
                        ntGamma[b1, e1, e3]}]], "ABORTED"],
  {True, True}];

(* --- D2. The EMITTED spelling of D. A dressed chain is a std::vector<DChainTok>: each fixed token is
   a DFac wrapped in dtfix(...), and dtr takes and returns a DFac, so a backward fixed token must read
   dtfix(dtr(DFac)). The reverse nesting, dtr(dtfix(...)), does not compile. --------------------- *)
NumTracer`Private`$ctCache = <||>;
check["dressed chain spells a backward fixed token dtfix(dtr(...))",
  CheckAbort[First @ First @ NumTracer`Private`compileDirac[
      {ntDeltaDirac[e1, e2], NumTracer`Private`ntDressedNum[{{1, {"ident"}}}, e3, e2],
       ntGamma[b1, e1, e3], ntVec[qq, b1]},
      <|b1 -> 0, e1 -> 100, e2 -> 101, e3 -> 102|>,
      <|qq -> <|"Base" -> 0, "Inv" -> None, "InvS" -> None|>|>, 15], "ABORTED"],
  "std::vector<DChainTok>{dtrslot(0), dtfix(dtr(dslash({{1.0,0}})))}"];

(* --- E. The same defect one level down, in a COLLECTED Dirac slot. A slot is an OPEN chain
   din->dout whose tokens are spliced into the surrounding loop IN CHAIN ORDER, so it needs exactly
   one open spinor label that is some head's `in` and no head's `out`. Here both open ends (e1, e2)
   are `in`, so there are two candidates; the old SelectFirst default picked one silently. The sum
   carries NON-numeric coefficients so it is genuinely collectible (diracSlotSumQ). ------------- *)
slot = Global`zA[Global`p] ntGamma[c1, e1, m1] ntGamma[c2, e2, m1] +
       Global`zB[Global`p] ntGamma[c1, e1, m2] ntGamma[c2, e2, m2];
check["anomalous collected slot is refused (slotorient)",
  CheckAbort[NumTracer`Private`diracSlotDecompose[slot], "ABORTED"],
  "ABORTED"];

If[fail, Print["orientation guard: FAILED"]; Exit[1]];
Print["orientation guard: all checks passed"];
