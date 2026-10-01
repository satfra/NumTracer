#!/usr/bin/env wolframscript
(* Codegen gate for a MULTI-TERM (Plus) Dirac projector of MIXED gamma parity: a hand-built DSL net
   (no FunKit / no TensorBases / no FORM / no dressings) reproducing the bug where the odd-trace
   filter in DSL.m judged parity on a DIAGRAM-GLOBAL gamma count and so discarded the whole diagram.

   The projector Plus below carries ALL-NUMERIC coefficients and lives in a single sector, so
   distributeQ is False (expandBridges leaves it EAGER for et::add) and collectibleDiracSumQ is
   False (that predicate is gated on dressedStructureSumQ, i.e. a NON-numeric coefficient). It is
   therefore neither distributed nor exempted — exactly the crack the bug fell through.

     proj = 1 + q3_slash                     branch gamma-counts 0 and 1  -> MIXED parity
     d1   = tr[ proj . q1_slash . ql_slash ]
          = tr[q1_slash ql_slash] + tr[q3_slash q1_slash ql_slash]
          = 4 (q1.ql)             + 0                                    (odd trace vanishes)

   The diagram-global gamma count is 0/1 + 1 + 1 = 3 -> ODD, so the pre-fix Select dropped the
   ENTIRE diagram, taking the non-vanishing EVEN branch with it: 0 diagrams, an identically-zero
   kernel, and no error raised. The correct verdict is per-branch: drop only when EVERY branch is
   odd. compare_multiterm_proj_num.cpp checks the emitted kernel against 4 (q1.ql) analytically,
   so a kernel that is merely non-empty but wrong still fails.

   The net depends on the loop angle through ql, so the surviving branch is genuinely exercised
   over the frame rather than folding to a constant. *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
                                   Global`q1, Global`q2, Global`q3, Global`ql];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p};

(* trailing '+' keeps the statement open across the line break (a leading '+' would end it) *)
proj = ntDeltaDirac[d1, d2] +
       ntGamma[c1, d1, d2] ntVec[Global`q3, c1];
net  = proj * ntGamma[a1, d2, d3] ntVec[Global`q1, a1] *
              ntGamma[a2, d3, d1] ntVec[Global`ql, a2];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];

(* Assertion A: the mixed-parity diagram must SURVIVE the odd-trace filter. Pre-fix this is 0. *)
If[Length[ntk[[1]]["Diagrams"]] < 1,
  Print["FAIL: the mixed-parity multi-term Dirac Plus was dropped by the odd-gamma filter ",
        "(0 diagrams). The filter must judge parity per BRANCH, not on the diagram-global ",
        "gamma count."]; Exit[1]];
Print["[gen] ok: mixed-parity multi-term projector survives the filter (",
      Length[ntk[[1]]["Diagrams"]], " diagram(s))"];

(* Assertion C (the other direction — the fix must not become a blanket exemption): a Plus whose
   branches are ALL odd once distributed is genuinely zero and must still be pruned.
     projZ = 1 + q1_slash q2_slash            branch gamma-counts 0 and 2
     dZ    = tr[ projZ . ql_slash ]           -> branch totals 1 and 3, BOTH odd -> identically 0
   (Its diagram-global count is also 3, so this case was pruned correctly by accident before; it is
   here to pin that the branch-aware rule keeps pruning it.) *)
projZ = ntDeltaDirac[e1, e2] +
        ntGamma[b1, e1, f1] ntVec[Global`q1, b1] ntGamma[b2, f1, e2] ntVec[Global`q2, b2];
netZ  = projZ * ntGamma[b3, e2, e1] ntVec[Global`ql, b3];
ntkZ  = NumTrace[netZ, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[Length[ntkZ[[1]]["Diagrams"]] =!= 0,
  Print["FAIL: an all-odd-branch diagram (identically zero) was NOT pruned (",
        Length[ntkZ[[1]]["Diagrams"]], " diagram(s)) — the odd-trace pruning has been ",
        "over-broadened into a blanket exemption."]; Exit[1]];
Print["[gen] ok: all-odd-branch diagram still pruned (0 diagrams)"];

(* Assertion D: POWERS must not blow up the parity fold. A scalar power like sp[q2,q3]^2 (ubiquitous
   in Gram determinants / projector normalisations) reaches diracParities as Power[b, n]. Folding it
   as `Times @@ ConstantArray[b, n]` recurses forever, because Times re-collapses the copies back to
   b^n — it hits $RecursionLimit and feeds garbage to Union, silently corrupting every verdict in the
   diagram. Same net as above, scaled by a squared scalar and a squared tensor self-contraction. *)
projP = ntDeltaDirac[g1, g2] + ntGamma[h1, g1, g2] ntVec[Global`q3, h1];
netP  = Global`p^2 * (1 + Global`p^2)^3 *
        projP * ntGamma[h2, g2, g3] ntVec[Global`q1, h2] *
                ntGamma[h3, g3, g1] ntVec[Global`ql, h3];
ntkP = Quiet@Check[NumTrace[netP, "Frame" -> frame, "Args" -> args, "Dressings" -> {}], $Failed];
If[ntkP === $Failed || Length[ntkP[[1]]["Diagrams"]] < 1,
  Print["FAIL: a diagram carrying Power[...] factors did not survive — the parity fold most likely ",
        "recursed on Power (Times re-collapses ConstantArray copies) and corrupted the verdict."];
  Exit[1]];
Print["[gen] ok: Power factors fold correctly (", Length[ntkP[[1]]["Diagrams"]], " diagram(s))"];

(* Assertion E: a COLLAPSED spinor loop must keep its tr(1) = 4. orderDiracFacs drops δ "connector"
   factors, so a loop built only from ntDeltaDirac arrives token-free; the runtime's split_loops
   discards empty segments, so without compileDirac's nEmptyLoops restoration the loop contributes
   nothing and the kernel comes out 4x too small per collapsed loop. (The dressed path already
   compensated for this; the bare path did not.) Emitted here and graded analytically in
   compare_multiterm_proj_num.cpp against 4 (tr[1]) * 3 (colour Nc) * (q1.ql). *)
netD = ntDeltaDirac[m1, m2] ntDeltaDirac[m2, m1] *
       ntSUNDeltaFund[3, n1, n2] ntSUNDeltaFund[3, n2, n1] *
       ntVec[Global`q1, r1] ntVec[Global`ql, r1];
ntkD = NumTrace[netD, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[Length[ntkD[[1]]["Diagrams"]] < 1,
  Print["FAIL: the collapsed-delta-loop diagram did not survive"]; Exit[1]];
MakeNTKernel[ntkD,
  FileNameJoin[{$gen, "gen_deltaloop_num.cpp"}],
  FileNameJoin[{$gen, "Deltaloop_num_kernel.hh"}],
  FileNameJoin[{$gen, "Deltaloop_num_kernels.hh"}],
  "Name" -> "Deltaloop_num_kernel", "Namespace" -> "deltaloop_num", "Dressings" -> {}];
Print["[gen] ok: collapsed delta-loop kernel emitted"];

(* Assertion F: a momentum-free pure-LORENTZ component must not be treated as a constant SU(N) one.
   "Constant" means "a constant SU(N) component" to Codegen's per-component dispatch, which hands it
   to compileColG — and compileColG only knows group heads. A closed metric loop g_{μν} g^{μν} = D is
   momentum-free but carries no colour, so it was emitted as
   `colFacG[ntMetric[v1,v2], <|v1 -> 0, ...|>]` raw into the generator .cpp (an instant clang error).
   Same defect as the collapsed-δ case above, one head-list short: EVERY non-SU(N) tensor head must
   appear in that FreeQ list, since it is not a momentum test. Graded analytically as 4 * (q1.ql). *)
netM = ntMetric[v1, v2] ntMetric[v2, v1] *
       ntVec[Global`q1, w1] ntVec[Global`ql, w1];
ntkM = NumTrace[netM, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[AnyTrue[ntkM[[1]]["Diagrams"][[1]]["Components"],
           #["Constant"] && ! FreeQ[#["Factors"], _ntMetric | _ntEpsilon] &],
  Print["FAIL: a pure-Lorentz component is flagged Constant and would be routed to compileColG, ",
        "which emits raw Mathematica into the generator."]; Exit[1]];
MakeNTKernel[ntkM,
  FileNameJoin[{$gen, "gen_metricloop_num.cpp"}],
  FileNameJoin[{$gen, "Metricloop_num_kernel.hh"}],
  FileNameJoin[{$gen, "Metricloop_num_kernels.hh"}],
  "Name" -> "Metricloop_num_kernel", "Namespace" -> "metricloop_num", "Dressings" -> {}];
Print["[gen] ok: momentum-free metric-loop kernel emitted"];

(* Assertion G: a CONSTANT colour/flavour component that is a SUM (a Plus of SU(N) heads).
   compileColG split a Times and expanded Power but had no Plus case, so the sum fell through to
   colFacG — which, alone among the per-head dispatchers in Codegen*.m, had neither a Plus branch nor
   a catch-all — and `colFacG[Plus[...], <|...|>]` was ToString'd raw into the generator .cpp.
   Found on the four-quark Fierz flavour structure delta*delta - 4*T*T; reproduced here WITHOUT
   FunKit so the cheap gate covers it.
   Such a sum is exempted by neither front-end predicate: sectorBridgeQ needs BOTH a colour and a
   Lorentz head, and dressedStructureSumQ needs a non-numeric coefficient. A pure-colour sum with
   numeric coefficients is therefore left eager, and the colour backend has no sum node at all
   (SUNNet is a flat product, no add(), no coefficient field) — so it must lower to a LIST of
   {net, scalar} branches that the generator's per-group colour sum adds up.
   Exact value: delta_{ij} delta_{ji} = Nc = 3, and T^a_{ij} T^a_{ji} = tr(T^a T^a) = (Nc^2-1)/2 = 4,
   so the sum is 3 + (-2)(4) = -5, times the Lorentz factor (q1.ql). Graded as -5 (q1.ql). *)
colSum = ntSUNDeltaFund[3, ci1, ci2] ntSUNDeltaFund[3, ci2, ci1] +
         (-2) ntSUNT[3, ca1, ci1, ci2] ntSUNT[3, ca1, ci2, ci1];
netG = colSum * ntVec[Global`q1, w2] ntVec[Global`ql, w2];
ntkG = NumTrace[netG, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
(* the colour sum must be ONE constant component, i.e. genuinely on the compileColG path *)
If[! AnyTrue[ntkG[[1]]["Diagrams"][[1]]["Components"],
             #["Constant"] && ! FreeQ[#["Factors"], _Plus] &],
  Print["FAIL: the pure-colour Plus is not a constant component, so this assertion no longer ",
        "covers the compileColG sum path. Components: ",
        ntkG[[1]]["Diagrams"][[1]]["Components"]]; Exit[1]];
MakeNTKernel[ntkG,
  FileNameJoin[{$gen, "gen_colsum_num.cpp"}],
  FileNameJoin[{$gen, "Colsum_num_kernel.hh"}],
  FileNameJoin[{$gen, "Colsum_num_kernels.hh"}],
  "Name" -> "Colsum_num_kernel", "Namespace" -> "colsum_num", "Dressings" -> {}];
Print["[gen] ok: constant colour-SUM kernel emitted"];

(* Assertion H: SU(N) FUNDAMENTAL Levi-Civita. NumTrace contracts epsilon PAIRS into Kronecker
   deltas (eps.eps = k! det(delta)) — a rewrite into existing primitives, like the gamma^0 one, so
   the C++ engine never sees an N-index object. Four CLOSED contractions, weighted by distinct
   primes so nothing cancels, each with an exactly hand-computable value:

     H1  SU(3), k=1, DIRECT closure:  eps[a,i1,i2] eps[a,j1,j2] d[i1,j1] d[i2,j2]
           = (d[i1,j1]d[i2,j2] - d[i1,j2]d[i2,j1]) d[i1,j1] d[i2,j2] = 3*3 - 3 =  6
     H2  SU(3), k=1, CROSSED closure: ... d[i1,j2] d[i2,j1]                       = 3 - 3*3 = -6
     H3  SU(2), k=0, DIRECT:  eps[F1,F3] eps[F2,F4] d[F1,F2] d[F3,F4]             = 2*2 - 2 =  2
     H4  SU(2), k=0, CROSSED: ... d[F1,F4] d[F3,F2]                               = 2 - 2*2 = -2

   The direct/crossed pairs differ ONLY by the antisymmetry sign, so a dropped Signature flips H2
   and H4 and moves the total from -10 to +22 — a sign error cannot hide as a small perturbation.
   Both ranks appear in one net, pinning that sunRankOf is respected per head.
     2*6 + 3*(-6) + 5*2 + 7*(-2) = 12 - 18 + 10 - 14 = -10, times the Lorentz factor (q1.ql). *)
epsC[i1_, i2_, j1_, j2_] := ntEpsFund[3, ea, i1, i2] ntEpsFund[3, ea, j1, j2];
epsF[f1_, f2_, f3_, f4_] := ntEpsFund[2, f1, f3] ntEpsFund[2, f2, f4];
h1 = epsC[c1, c2, c3, c4] ntSUNDeltaFund[3, c1, c3] ntSUNDeltaFund[3, c2, c4];
h2 = epsC[e1, e2, e3, e4] ntSUNDeltaFund[3, e1, e4] ntSUNDeltaFund[3, e2, e3];
h3 = epsF[g1, g2, g3, g4] ntSUNDeltaFund[2, g1, g2] ntSUNDeltaFund[2, g3, g4];
h4 = epsF[k1, k2, k3, k4] ntSUNDeltaFund[2, k1, k4] ntSUNDeltaFund[2, k3, k2];
netH = (2 h1 + 3 h2 + 5 h3 + 7 h4) ntVec[Global`q1, w3] ntVec[Global`ql, w3];
ntkH = NumTrace[netH, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
(* H-a: the rewrite actually ran — no epsilon may survive into the analysed kernel. *)
If[! FreeQ[ntkH, _ntEpsFund],
  Print["FAIL: ntEpsFund survived NumTrace; expandFundEps did not run (or missed a context). ",
        "It would then reach colFacG, which has no rule for it."]; Exit[1]];
Print["[gen] ok: epsilon pairs contracted away (no ntEpsFund survives)"];
MakeNTKernel[ntkH,
  FileNameJoin[{$gen, "gen_epsfund_num.cpp"}],
  FileNameJoin[{$gen, "Epsfund_num_kernel.hh"}],
  FileNameJoin[{$gen, "Epsfund_num_kernels.hh"}],
  "Name" -> "Epsfund_num_kernel", "Namespace" -> "epsfund_num", "Dressings" -> {}];
Print["[gen] ok: fundamental-epsilon kernel emitted"];

(* Assertion H-c/H-d: the loud refusals. These Abort[], so CheckAbort is required — Quiet@Check
   (used elsewhere in this file) does NOT catch an Abort. *)
If[CheckAbort[NumTrace[ntEpsFund[3, z1, z2, z3] ntSUNDeltaFund[3, z1, z2] ntVec[Global`q1, z9] ntVec[Global`ql, z9],
     "Frame" -> frame, "Args" -> args, "Dressings" -> {}]; False, True] =!= True,
  Print["FAIL: an UNPAIRED fundamental epsilon was accepted. A lone epsilon is not an SU(N) ",
        "invariant and has no engine representation; it must abort (NumTrace::epsodd)."]; Exit[1]];
Print["[gen] ok: unpaired epsilon refused"];
(* distinct labels on purpose: identical ones would auto-collapse to Power[eps,2], which is the
   self-contraction path, not the arity path we mean to exercise here. *)
If[CheckAbort[NumTrace[ntEpsFund[3, y1, y2] ntEpsFund[3, y3, y4] ntVec[Global`q1, y9] ntVec[Global`ql, y9],
     "Frame" -> frame, "Args" -> args, "Dressings" -> {}]; False, True] =!= True,
  Print["FAIL: a rank-3 epsilon with only 2 indices was accepted (NumTracer::epsrank). The ",
        "fundamental epsilon of SU(N) carries exactly N indices."]; Exit[1]];
Print["[gen] ok: wrong-arity epsilon refused"];

(* Assertion J: an epsilon pair STRADDLING an eager Plus, with the two ranks INTERLEAVED. This is
   the shape of the diquark vertex of the four-quark Fierz bases:

     eps3[a,j1,j2] eps2[f1,f2] * ( 2 eps3[a,j3,j4] eps2[f3,f4] D1 - ... )

   Neither partner of either pair is reachable at one multiplicative level: the outer epsilons sit
   OUTSIDE the sum, their partners inside it. Each summand on its own holds exactly one rank-3 and
   one rank-2 epsilon, so a rank-blind pairing joins THOSE two — building a 2x2 determinant out of a
   3-index object and silently dropping a colour leg. That is not a small error: the two summands
   then expose DIFFERENT free indices and the run dies far downstream in checkLabels with
   NumTrace::plusfree, pointing at the sum rather than at the epsilon.

   The fix is two-part and this one assertion pins both, since reverting EITHER half breaks it:
     * pairing restricted to equal rank  -> without it, the mispairing above (wrong number);
     * distribution into the sum         -> without it, all four survive and epsodd aborts.

   Oracle: each summand factorises into the SU(3) and SU(2) contractions already pinned by H,
     summand 1 = 2 * (direct colour = 6) * (direct flavour = 2)  =  24
     summand 2 = 3 * (crossed colour = -6) * (crossed flavour = -2) = 36
   so the total is 60 * (q1.ql). The crossed summand carries TWO sign flips, so it is graded by the
   product: a single dropped Signature moves it to -36 and the total to -12, not to a near miss. *)
netJ = ntEpsFund[3, ja, j1, j2] ntEpsFund[2, jf1, jf2] *
  (2 ntEpsFund[3, ja, j3, j4] ntEpsFund[2, jf3, jf4] *
       ntSUNDeltaFund[3, j1, j3] ntSUNDeltaFund[3, j2, j4] *
       ntSUNDeltaFund[2, jf1, jf3] ntSUNDeltaFund[2, jf2, jf4] +
   3 ntEpsFund[3, ja, j5, j6] ntEpsFund[2, jf5, jf6] *
       ntSUNDeltaFund[3, j1, j6] ntSUNDeltaFund[3, j2, j5] *
       ntSUNDeltaFund[2, jf1, jf6] ntSUNDeltaFund[2, jf2, jf5]) *
  ntVec[Global`q1, w5] ntVec[Global`ql, w5];
ntkJ = NumTrace[netJ, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[! FreeQ[ntkJ, _ntEpsFund],
  Print["FAIL: an ntEpsFund survived a pair straddling an eager Plus — the sum was not ",
        "distributed, so the partners never met at one multiplicative level."]; Exit[1]];
Print["[gen] ok: epsilon pair straddling an eager Plus contracted away"];
MakeNTKernel[ntkJ,
  FileNameJoin[{$gen, "gen_epsstraddle_num.cpp"}],
  FileNameJoin[{$gen, "Epsstraddle_num_kernel.hh"}],
  FileNameJoin[{$gen, "Epsstraddle_num_kernels.hh"}],
  "Name" -> "Epsstraddle_num_kernel", "Namespace" -> "epsstraddle_num", "Dressings" -> {}];
Print["[gen] ok: straddling-epsilon kernel emitted"];

(* J-b: two epsilons of DIFFERENT rank and no same-rank partner for either. They must NOT be paired
   with each other (that is the silent-garbage path); they are left alone and reported by epsodd. *)
If[CheckAbort[NumTrace[ntEpsFund[3, m1, m2, m3] ntEpsFund[2, m4, m5] *
     ntSUNDeltaFund[3, m1, m2] ntSUNDeltaFund[2, m4, m5] ntVec[Global`q1, m9] ntVec[Global`ql, m9],
     "Frame" -> frame, "Args" -> args, "Dressings" -> {}]; False, True] =!= True,
  Print["FAIL: a rank-3 epsilon was contracted against a rank-2 one. epsPairExpand would size the ",
        "determinant by one partner and silently DROP the other's surplus indices."]; Exit[1]];
Print["[gen] ok: cross-rank epsilon pairing refused"];

(* J-c: four epsilons of the SAME rank, none sharing an index. Equal rank does not identify the
   index SPACE (at Nc == Nf a colour and a flavour epsilon are indistinguishable), and with nothing
   shared the partnering is not determined, so NumTracer must refuse rather than guess. The FORCED
   case — exactly two of a rank, as in J above — stays legal and is covered there. *)
If[CheckAbort[NumTrace[ntEpsFund[3, u1, u2, u3] ntEpsFund[3, u4, u5, u6] *
     ntEpsFund[3, u7, u8, u9] ntEpsFund[3, v1, v2, v3] ntVec[Global`q1, v9] ntVec[Global`ql, v9],
     "Frame" -> frame, "Args" -> args, "Dressings" -> {}]; False, True] =!= True,
  Print["FAIL: an ambiguous zero-overlap epsilon pairing was guessed at (NumTrace::epsambig). ",
        "Pairing across index spaces contracts colour against flavour and returns a wrong number."];
  Exit[1]];
Print["[gen] ok: ambiguous zero-overlap epsilon pairing refused"];

(* Assertion K: a CLOSED colour loop inside a colour-ENTANGLED Plus, i.e. a group head that reaches
   splitColourGroupsInv as an integer POWER rather than as a bare factor.

   splitColourGroupsInv distributes the entangled Pluses and then splits each branch into a colour
   product and a Lorentz/Dirac remainder with a LEVEL-1 Cases/DeleteCases, which matches only what is
   a BARE factor. A group head raised to a power is not, so it stays in the REMAINDER and is handed
   to the Lorentz-only builderInv, which has no rule for it and CForm's it into the generated C++.

   TWO entangled sums are required to reproduce, and this is the whole subtlety. With one sum the
   branch's factors reach the split as a LIST (Join of the term's factors and the un-expanded ones),
   so two identical deltas stay two separate elements and the head pattern matches both. It is
   Expand[Times @@ needExpand] across SEVERAL sums that multiplies them inside one term, where Times
   collapses them to deltaFund[3,ck1,ck2]^2 — a closed colour loop, value Nc.

   Both sums must also carry the colour AND the Dirac structure: with colour in a separate factor the
   two land in disconnected components (they share no label) and colour goes down the constant /
   compileColGSum path, which has always expanded such powers.

   Oracle — every one of the 4 branches has colour deltaFund^2 = Nc = 3; the Dirac side closes
   through ntDeltaDirac[dk2,dk1]:
     (2,3): 6  * tr(1)               = 24        (7,3),(2,5): odd gamma count -> 0
     (7,5): 35 * tr(q1slash qlslash) = 140*(q1.ql)
   so the kernel is 3*(24 + 140*(q1.ql)) = 72 + 420*(q1.ql). Constant and momentum-dependent parts
   are graded together, so dropping either surviving branch is visible. *)
netK = (2 ntDeltaDirac[dk1, dkx] ntSUNDeltaFund[3, ck1, ck2] +
        7 ntGamma[mk1, dk1, dkx] ntVec[Global`q1, mk1] ntSUNDeltaFund[3, ck1, ck2]) *
       (3 ntDeltaDirac[dkx, dk2] ntSUNDeltaFund[3, ck1, ck2] +
        5 ntGamma[mk2, dkx, dk2] ntVec[Global`ql, mk2] ntSUNDeltaFund[3, ck1, ck2]) *
       ntDeltaDirac[dk2, dk1];
ntkK = NumTrace[netK, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
MakeNTKernel[ntkK,
  FileNameJoin[{$gen, "gen_colpow_num.cpp"}],
  FileNameJoin[{$gen, "Colpow_num_kernel.hh"}],
  FileNameJoin[{$gen, "Colpow_num_kernels.hh"}],
  "Name" -> "Colpow_num_kernel", "Namespace" -> "colpow_num", "Dressings" -> {}];
Print["[gen] ok: closed-colour-loop (deltaFund^2) in an entangled Plus emitted"];

(* Assertion I: the ADJOINT Levi-Civita. !! SU(2) ONLY !! At rank 2 — and ONLY at rank 2 — the
   adjoint epsilon coincides with the structure constant, eps^abc = f^abc exactly (coefficient +1),
   so FromFunKit rewrites it to ntSUNf[2,...]. This is NOT an SU(N) statement: for SU(3),
   f^123 = 1 but f^147 = 1/2, so f is not an epsilon at all there, and rank != 2 is REFUSED.

   The +1 is verified independently of NumTracer (Pauli matrices, T^a = sigma^a/2,
   f = -2i tr([T^a,T^b]T^c), matching sun_net.hpp:225,233-247): max|f - eps| = 0.

   Graded with a contraction LINEAR in the coefficient — eps^abc f^abc = f^abc f^abc = 3! = 6 for
   SU(2). Every contraction QUADRATIC in eps is invariant under c -> -c, so it would pass with the
   sign inverted; only a linear one pins the sign. 6 * (q1.ql).
   (FunKit itself is not loaded here: FromFunKit matches these heads by SymbolName, so plain Global
   symbols exercise the same rewrite.) *)
netI = FromFunKit[Global`epsAdjFlav[ia, ib, ic] Global`fFlav[ia, ib, ic]] *
       ntVec[Global`q1, w4] ntVec[Global`ql, w4];
If[! FreeQ[netI, Global`epsAdjFlav],
  Print["FAIL: epsAdjFlav was not rewritten by FromFunKit."]; Exit[1]];
ntkI = NumTrace[netI, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
MakeNTKernel[ntkI,
  FileNameJoin[{$gen, "gen_epsadj_num.cpp"}],
  FileNameJoin[{$gen, "Epsadj_num_kernel.hh"}],
  FileNameJoin[{$gen, "Epsadj_num_kernels.hh"}],
  "Name" -> "Epsadj_num_kernel", "Namespace" -> "epsadj_num", "Dressings" -> {}];
Print["[gen] ok: SU(2) adjoint-epsilon kernel emitted (eps = f, SU(2) ONLY)"];

(* Assertion M: the refusals around the adjoint epsilon and the untranslated-head class.
   M1 — adjoint epsilon at rank != 2. Global`Nc = 3 here, so epsAdjCol must abort: the SU(3)
   adjoint epsilon carries 8 indices, not 3, and ShowFormTracerDefinitions[]'s generic
   epsAdjCol[a,b,c] display is exactly what would mislead someone into writing it. *)
If[CheckAbort[FromFunKit[Global`epsAdjCol[ja, jb, jc] Global`FCol[ja, jb, jc]]; False, True] =!= True,
  Print["FAIL: an adjoint epsilon at SU(3) was accepted. eps = f holds ONLY for SU(2); at SU(3) ",
        "the adjoint epsilon has 8 indices."]; Exit[1]];
Print["[gen] ok: SU(3) adjoint epsilon refused (eps = f is SU(2)-only)"];
(* M2 — the ROOT-CLASS guard: a declared FunKit head with no map entry must abort, not become a
   silent scalar coefficient. epsLorentz is the 3D SPATIAL epsilon; ntEpsilon is 4D and hard-wired
   to four labels, so mapping one onto the other would be a silent dimension error. This is the
   guard that would have caught epsFundCol immediately. *)
If[CheckAbort[FromFunKit[Global`epsLorentz[la, lb, lc] ntVec[Global`q1, la]]; False, True] =!= True,
  Print["FAIL: epsLorentz (a declared FunKit head with no map entry) was accepted. An untranslated ",
        "head becomes a SCALAR COEFFICIENT with invisible indices — the exact bug this guards."];
  Exit[1]];
Print["[gen] ok: untranslated FunKit head refused (root-class guard)"];

Print["[gen] Numeric backend ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_multiterm_proj_num.cpp"}],
  FileNameJoin[{$gen, "Multiterm_proj_num_kernel.hh"}],
  FileNameJoin[{$gen, "Multiterm_proj_num_kernels.hh"}],
  "Name" -> "Multiterm_proj_num_kernel", "Namespace" -> "multiterm_proj_num",
  "Dressings" -> {}];

Print["Multi-term Dirac projector codegen kernels generated."];
