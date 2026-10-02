(* ::Package:: *)

(* DSL analysis: classify the heads, split each term into independent contraction
   components, and allocate the env-id layout the code generator consumes.

   Loaded inside NumTracer`Private` by NumTracer.m — public symbols (NumTrace, the
   nt* heads) already exist in the NumTracer` context. *)

(* ---- environment flags ------------------------------------------------------- *)

(* The ONE truth test for every NT_* boolean env flag (here because DSL.m loads first and the
   Codegen*.m files share this private context). Only 1/true/yes/on read as ON, so `VAR=0` is OFF,
   and an unset variable (Environment[] = $Failed, e.g. after SetEnvironment["VAR" -> None]) is OFF. *)
ntEnvFlag[name_String] :=
  With[{v = Environment[name]},
    StringQ[v] && MemberQ[{"1", "true", "yes", "on"}, ToLowerCase[StringTrim[v]]]];

(* A positive-integer environment variable, or 0 for "unset / unusable" (0 is also the manifest's
   spelling of "no cap"). Parsed as digits, never evaluated. *)
ntEnvPosInt[name_String] :=
  With[{v = Environment[name]},
    If[StringQ[v] && StringMatchQ[StringTrim[v], DigitCharacter ..], FromDigits[StringTrim[v]], 0]];

(* ---- [prof] accumulators -----------------------------------------------------
   ntProfTimed[key, expr] evaluates expr. While $ntProfOn is True it also adds the wall time to
   $ntProf[key] = {calls, inclusive s, exclusive s}: a key already on the stack (recursion) is timed at
   its outermost call only, and exclusive time excludes nested timed keys. Callers Block $ntProfOn and
   $ntProf around a stage and print the result with ntProfReport. Off, it costs ~0.3 µs per call; on,
   ~9 µs, which is visible on stages with ~10^5 timed calls (compileDirac). *)
$ntProfOn = False;
$ntProf = <||>;
$ntProfStack = {};
$ntProfChild = 0.;
SetAttributes[{ntProfTimed, ntProfRun}, HoldRest];
ntProfTimed[key_, expr_] := If[$ntProfOn, ntProfRun[key, expr], expr];
ntProfRun[key_, expr_] :=
  If[MemberQ[$ntProfStack, key],
    expr,
    Module[{t, res, child},
      Block[{$ntProfStack = Append[$ntProfStack, key], $ntProfChild = 0.},
        {t, res} = AbsoluteTiming[expr];
        child = $ntProfChild];
      $ntProfChild += t;
      $ntProf[key] = Lookup[$ntProf, key, {0, 0., 0.}] + {1, t, t - child};
      res]];

ntProfReport[prefix_String] :=
  KeyValueMap[
    ntLog[prefix, #1, ": ", #2[[1]], " calls, ", #2[[2]], " s incl, ", #2[[3]], " s excl"] &,
    ReverseSortBy[$ntProf, #[[2]] &]];

(* ---- head registry ------------------------------------------------------------

   ONE record per tensor head (a factor that takes part in the contraction, vs. a scalar
   coefficient). Every head classification below is derived from this list at load time, so a new
   head is registered here and nowhere else. Keys:
     "Form"      (held) the head with named argument patterns; the accessors are defined on it
     "Sector"    "Lorentz" | "Dirac" | "Adjoint" | "Fundamental". Lorentz, Dirac and the SU(N) groups
                 contract in disjoint index spaces. The SU(N) heads (Adjoint, Fundamental) carry their
                 rank N as the FIRST argument, so several groups coexist in one network
     "Labels"    (held, in Form's names) the contraction index labels                -> labelsOf
     "Spinor"    (held) the subset of Labels in the spinor index space (default {})  -> spinorLabelsHead
     "Momentum"  the first argument is the momentum the head carries (the source of runtime Var
                 leaves); the value is the head's sign under q -> -q                 -> momentumOf
     "Inv"       the momentum needs a 1/q^2 env slot (a projector with a full denominator)
     "InvS"      the momentum needs a SPATIAL 1/|q⃗|^2 env slot (finite-T electric/magnetic projectors)
     "Tensor"    the tensorQ pattern, where narrower than the bare head
     "Rewritten" NumTrace expands the head away (expandFundEps), so it is not a codegen colour head
   The order is the order of the derived Alternatives. *)
$ntHeads = {
  <|"Form" :> ntMetric[mu_, nu_], "Sector" -> "Lorentz", "Labels" :> {mu, nu}|>,
  (* ntVec[q, i_Integer] is the scalar component q_i (0-based, 0 = temporal), resolved by the frame
     like ntSP, so it is NOT a tensor. *)
  <|"Form" :> ntVec[q_, mu_], "Sector" -> "Lorentz", "Labels" :> {mu}, "Momentum" -> -1,
    "Tensor" -> ntVec[_, Except[_Integer]]|>,
  <|"Form" :> ntTransProj[q_, mu_, nu_], "Sector" -> "Lorentz", "Labels" :> {mu, nu}, "Momentum" -> 1,
    "Inv" -> True|>,
  <|"Form" :> ntLongProj[q_, mu_, nu_], "Sector" -> "Lorentz", "Labels" :> {mu, nu}, "Momentum" -> 1,
    "Inv" -> True|>,
  (* P_E = P_T − P_M uses both 1/q² and 1/|q⃗|² *)
  <|"Form" :> ntElectricProj[q_, mu_, nu_], "Sector" -> "Lorentz", "Labels" :> {mu, nu}, "Momentum" -> 1,
    "Inv" -> True, "InvS" -> True|>,
  <|"Form" :> ntMagneticProj[q_, mu_, nu_], "Sector" -> "Lorentz", "Labels" :> {mu, nu}, "Momentum" -> 1,
    "InvS" -> True|>,
  <|"Form" :> ntSUNf[n_, a_, b_, c_], "Sector" -> "Adjoint", "Labels" :> {a, b, c}|>,
  <|"Form" :> ntSUNDeltaAdj[n_, a_, b_], "Sector" -> "Adjoint", "Labels" :> {a, b}|>,
  <|"Form" :> ntGamma[mu_, din_, dout_], "Sector" -> "Dirac", "Labels" :> {mu, din, dout},
    "Spinor" :> {din, dout}|>,
  <|"Form" :> ntGamma5[din_, dout_], "Sector" -> "Dirac", "Labels" :> {din, dout}, "Spinor" :> {din, dout}|>,
  <|"Form" :> ntC[din_, dout_], "Sector" -> "Dirac", "Labels" :> {din, dout}, "Spinor" :> {din, dout}|>,
  (* ntSigma's OPEN Lorentz labels are its FREE legs (gluon ids); slashed legs carry a momentum, not
     an open id. The spinor axes din,dout are always open. *)
  <|"Form" :> ntSigma[legA_, legB_, din_, dout_], "Sector" -> "Dirac",
    "Labels" :> Join[Cases[{legA, legB}, {"free", mu_} :> mu], {din, dout}], "Spinor" :> {din, dout}|>,
  <|"Form" :> ntDeltaDirac[din_, dout_], "Sector" -> "Dirac", "Labels" :> {din, dout},
    "Spinor" :> {din, dout}|>,
  <|"Form" :> ntSUNT[n_, a_, i_, j_], "Sector" -> "Fundamental", "Labels" :> {a, i, j}|>,
  <|"Form" :> ntSUNDeltaFund[n_, i_, j_], "Sector" -> "Fundamental", "Labels" :> {i, j}|>,
  (* the Lorentz Levi-Civita ε_{μνρσ} *)
  <|"Form" :> ntEpsilon[a_, b_, c_, d_], "Sector" -> "Lorentz", "Labels" :> {a, b, c, d}|>,
  (* per-component-dressed δ (see their ::usage): classified exactly like the plain δ; the dressing
     is folded numerically by sun_value_dressed *)
  <|"Form" :> ntSUNDiagFund[n_, i_, j_, spec_], "Sector" -> "Fundamental", "Labels" :> {i, j}|>,
  <|"Form" :> ntSUNDiagAdj[n_, a_, b_, spec_], "Sector" -> "Adjoint", "Labels" :> {a, b}|>,
  (* SU(N) fundamental Levi-Civita: N indices, all of them contraction labels *)
  <|"Form" :> ntEpsFund[n_, idx___], "Sector" -> "Fundamental", "Labels" :> {idx}, "Rewritten" -> True|>,
  (* ntDressedNum[options, din, dout]: a dressed propagator numerator kept EAGER (symbolic dressing
     collection). `options` is a tensor-FREE list of {coeffExpr, spec} (spec = {"ident"} |
     {"slash",vlc}), so it exposes only its spinor labels; the dressings fold in C++ (one DPoly trace)
     instead of distributing the numerator into 2^D diagrams. *)
  <|"Form" :> ntDressedNum[opts_, din_, dout_], "Sector" -> "Dirac", "Labels" :> {din, dout},
    "Spinor" :> {din, dout}|>,
  (* ntDiracSlot[options, din, dout, legs]: a collected Dirac slot, a coefficient-weighted sum of
     Dirac structures sharing the spinor pair (din,dout) AND the same set of open Lorentz legs `legs`
     (k=0: a propagator numerator, k>=1: a vertex). `options` is a list of {coeffExpr,
     structureProduct} whose free Lorentz ids are exactly `legs`, so the surrounding net closes them
     for every structure: it exposes din,dout AND every open leg as contraction labels. *)
  <|"Form" :> ntDiracSlot[opts_, din_, dout_, legs_], "Sector" -> "Dirac",
    "Labels" :> Join[legs, {din, dout}], "Spinor" :> {din, dout}|>
};

ntHeadSym[rec_] := Replace[Extract[rec, Key["Form"], Hold], Hold[h_Symbol[___]] :> h];
ntTensorPatOf[rec_] := Lookup[rec, "Tensor", Blank[ntHeadSym[rec]]];
ntSectorHeads[sectors_List] := Select[$ntHeads, MemberQ[sectors, #["Sector"]] &];
ntSectorPat[sectors_List] := Alternatives @@ (ntTensorPatOf /@ ntSectorHeads[sectors]);

(* accessor[Form] := rec[key], with the held key value as the body *)
ntDefineAccessor[accessor_Symbol, rec_, key_String] :=
  Replace[{Extract[rec, Key["Form"], Hold], Extract[rec, Key[key], Hold]},
    {Hold[lhs_], Hold[rhs_]} :> (accessor[lhs] := rhs)];

$ntTensorPat      = Alternatives @@ (ntTensorPatOf /@ $ntHeads);
$ntAdjointPat     = ntSectorPat[{"Adjoint"}];
$ntFundamentalPat = ntSectorPat[{"Fundamental"}];
$ntSUNHeadPat     = ntSectorPat[{"Adjoint", "Fundamental"}];
(* by bare head, so a component carrying an integer component ntVec[q, i] is not constant either *)
$ntNonSUNHeadPat  = Alternatives @@ (Blank @* ntHeadSym /@ ntSectorHeads[{"Lorentz", "Dirac"}]);

(* ---- head classification ---------------------------------------------------- *)

With[{p = $ntTensorPat}, tensorQ[p] = True];
tensorQ[_] = False;

With[{p = $ntAdjointPat}, adjointSUNQ[p] = True];          (* group-adjoint heads (bridge with Lorentz) *)
adjointSUNQ[_] = False;
With[{p = $ntFundamentalPat}, fundamentalSUNQ[p] = True];  (* group-fundamental heads (quark-line) *)
fundamentalSUNQ[_] = False;

(* The SU(N) rank a group head builds against: its leading argument. *)
sunRankOf[h_] := First[h];

(* labelsOf: the contraction index labels of a tensor factor (momentum arg dropped).
   spinorLabelsHead: its spinor (Dirac) axis labels, which get a disjoint id range so the engine never
   contracts a spinor axis against a Lorentz/colour axis sharing an id.
   momentumOf: the momentum it carries, or None.
   needsInvQ / needsInvSQ: whether that momentum needs a 1/q^2 / spatial 1/|q⃗|^2 env slot. *)
Scan[Function[rec,
    ntDefineAccessor[labelsOf, rec, "Labels"];
    If[KeyExistsQ[rec, "Spinor"], ntDefineAccessor[spinorLabelsHead, rec, "Spinor"]];
    If[KeyExistsQ[rec, "Momentum"],
      Replace[Extract[rec, Key["Form"], Hold],
        Hold[form : _[Verbatim[Pattern][mom_Symbol, _], ___]] :> (momentumOf[form] := mom)]];
    With[{h = ntHeadSym[rec]},
      If[TrueQ[rec["Inv"]], needsInvQ[h[__]] = True];
      If[TrueQ[rec["InvS"]], needsInvSQ[h[__]] = True]]],
  $ntHeads];
spinorLabelsHead[_] := {};
momentumOf[_]       := None;
needsInvQ[_]        = False;
needsInvSQ[_]       = False;

(* All spinor labels anywhere under a (sub)expression. *)
allSpinorLabels[e_] := DeleteDuplicates @ Flatten @ Cases[e, h_?tensorQ :> spinorLabelsHead[h], {0, Infinity}];

(* ---- self-trace normalization ----------------------------------------------- *)

(* A tensor with a repeated index is a self-trace (e.g. P^mu_mu). The engine
   contracts pairwise BETWEEN tensors and never self-contracts one, so we relabel
   the second occurrence and insert the matching identity (Lorentz metric for a
   Lorentz index, adjoint delta for a colour index): P^mu_mu = P^{mu nu} d_{mu nu}.
   This keeps every index appearing on two distinct tensors, as the engine needs. *)
(* A SUM vertex (tensorQ[Plus] is False) passes through untouched, so a repeated index INSIDE a
   summand is never split. No flow produces one, and labelCensus's Plus branch would flag it; fixing
   it here would need every summand's free-index set kept aligned, which this local rewrite cannot see. *)
splitSelfTraces[factors_List] := Module[{res = {}, conns = {}},
  Function[f, If[! tensorQ[f], AppendTo[res, f],
    (* The connecting identity reuses the SAME group rank N as the head it closes: an
       adjoint group self-trace closes with ntSUNDeltaAdj[N,..], a fundamental one with
       ntSUNDeltaFund[N,..], a Lorentz/Dirac one with the metric. *)
    Module[{dups, relabeled = f, conn = Which[
        adjointSUNQ[f],     With[{n = sunRankOf[f]}, ntSUNDeltaAdj[n, ##] &],
        fundamentalSUNQ[f], With[{n = sunRankOf[f]}, ntSUNDeltaFund[n, ##] &],
        True,               ntMetric]},                          (* Lorentz/Dirac (unchanged) *)
      dups = Cases[Tally[labelsOf[f]], {l_, c_} /; c >= 2 :> l];
      Do[With[{fresh = Unique["st"]},
           relabeled = ReplacePart[relabeled, Last[Position[relabeled, l, {1}]] -> fresh];  (* relabel 2nd occurrence *)
           AppendTo[conns, conn[l, fresh]]], {l, dups}];
      AppendTo[res, relabeled]]]] /@ factors;
  Join[res, conns]
];

(* ---- free indices & scalar test (work through Plus/Times, for eager summation) ---- *)

(* Whether a (sub)expression carries no tensor head — a pure scalar coefficient. *)
With[{p = $ntTensorPat}, scalarQ[e_] := FreeQ[e, p]];

(* The free (uncontracted) index labels of a tensor (sub)expression. A product sums
   indices that appear twice (free = appear once); a sum's summands share free indices
   (a vertex's legs); a bare head exposes all its labels. Used to group components and
   to align eager add(...) operands. *)
freeIdx[e_] := Which[
  tensorQ[e],        labelsOf[e],
  Head[e] === Plus,  freeIdx[First[List @@ e]],
  (* A tensor^n is n copies sharing the SAME labels — a closed self-contraction (see compileLorentz),
     so it exposes NO free index. Spelled out rather than left to the True branch below, which
     would return {} for the wrong reason and hide a malformed Power. *)
  Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 2 && ! scalarQ[e], {},
  (* count == 1 is free; count == 2 is contracted (checkLabels guarantees no label occurs more
     often). Not OddQ, which would misclassify a 3x/4x label if that guard were bypassed. *)
  Head[e] === Times, Cases[Tally[Flatten[freeIdx /@ (List @@ e)]], {l_, c_} /; c == 1 :> l],
  True,              {}];

(* All index labels anywhere under a (sub)expression (free + internally summed) — the
   set that needs an axis id. *)
allLabels[e_] := DeleteDuplicates @ Flatten @ Cases[e, h_?tensorQ :> labelsOf[h], {0, Infinity}];

(* ---- connected components over shared free indices -------------------------- *)

(* Group tensor factors (heads OR Plus-vertices) into independent sub-networks: two
   factors are connected if their free indices intersect. Each component contracts to
   a scalar (colour stays a separate constant); keeps every contraction small. *)
connectedComponents[factors_List] := Module[{byLabel, edges},
  (* label -> the factors carrying it; every pair of those is an edge. Union keeps the edge list in
     the lexicographic order Subsets gives, so the Graph (and its component order) is unchanged. *)
  byLabel = GroupBy[
    Join @@ MapIndexed[Function[{ls, i}, {#, First[i]} & /@ ls], freeIdx /@ factors],
    First -> Last];
  edges = Union @@ (Subsets[Union[#], {2}] & /@ Values[byLabel]);
  factors[[#]] & /@ ConnectedComponents[Graph[Range[Length[factors]], UndirectedEdge @@@ edges]]
];

(* Greedy contraction order: keep each successive factor sharing a free index with the
   running set, so contract_all's intermediates stay low-rank (never outer-product the
   whole thing). *)
orderFactors[fs_List] := Module[{fi = freeIdx /@ fs, byLabel, touched, rem, out, take},
  If[fs === {}, Return[{}]];
  (* touched[[i]]: factor i shares a free index with the factors taken so far *)
  byLabel = GroupBy[Join @@ MapIndexed[Function[{ls, i}, {#, First[i]} & /@ ls], fi], First -> Last];
  touched = ConstantArray[False, Length[fs]];
  take = Function[k, AppendTo[out, k]; Scan[(touched[[#]] = True) &, Join @@ Lookup[byLabel, fi[[k]]]]];
  out = {}; take[1]; rem = Range[2, Length[fs]];
  While[rem =!= {},
    With[{pick = SelectFirst[rem, touched[[#]] &, First[rem]]},
      take[pick]; rem = DeleteCases[rem, pick, {1}, 1]]];
  fs[[out]]];

(* ---- sector-bridge expansion (keep colour and Lorentz contractions separate) ---- *)

(* Eager summation (an add(...) net) is a win only when a structure-sum lives in ONE sector (e.g. a
   3-gluon vertex: Lorentz structures times one colour factor). A SECTOR-BRIDGING sum pairs each
   colour structure with a Lorentz one (e.g. a 4-gluon vertex's f.f ⊗ metric.metric); add(...)-ing it
   fuses colour and Lorentz axes into one tensor whose entry count explodes. So DISTRIBUTE it: each
   term is sector-disjoint again, at a cost linear in the (few) structures. *)
(* The Lorentz side excludes ntEpsilon: an ε-carrying colour sum stays eager. *)
With[{adj = $ntAdjointPat, lor = DeleteCases[ntSectorPat[{"Lorentz"}], Verbatim[_ntEpsilon]]},
  sectorBridgeQ[p_] := (! FreeQ[p, adj]) && (! FreeQ[p, lor])];

(* The scalar (non-tensor) coefficient of a single summand. *)
scalarCoeffOf[t_] := Times @@ Select[If[Head[t] === Times, List @@ t, {t}], scalarQ];
(* A TENSOR-structure sum whose per-summand scalar coefficients are not all numeric must be
   DISTRIBUTED too: the eager add(...) net can only carry NUMERIC per-structure scalars (each
   summand is scaled by a numeric literal). A fermion propagator numerator
   Mq*deltaDirac + (.. dressings ..)*(gamma.vec) carries runtime dressing coefficients, so
   each Dirac structure must become its own diagram with that dressing as a scalar coeff.
   (Pure-scalar dressing sums — no tensor head — are left intact as coefficients.) *)
dressedStructureSumQ[p_Plus] := (! scalarQ[p]) && AnyTrue[List @@ p, ! NumericQ[scalarCoeffOf[#]] &];
(* Symbolic dressing collection: when $ntDressCollect is True a collectible DIRAC dressed structure
   sum (e.g. a propagator numerator Mq·δ + Z(p)·γ·p) is kept eager and rewritten to one ntDressedNum /
   ntDiracSlot token (rewriteDressedNums), so the diagram folds to ONE trace carrying its dressings
   instead of exploding into 2^D diagrams. Sector-bridging and non-Dirac dressed sums still distribute.
   NumTrace and FromFunKit set it from their "DressingCollection" option (default True). *)
$ntDressCollect = False;
(* The OPEN (free) spinor labels of a summand: spinor indices appearing an odd number of times. *)
openSpinorOf[t_] := Cases[Tally[Flatten[Cases[t, h_?tensorQ :> spinorLabelsHead[h], {0, Infinity}]]],
  {l_, c_} /; OddQ[c] :> l];
(* A propagator-numerator structure sum: every term carries Dirac structure and the SAME pair of open
   spinor indices (so it connects one spinor-in to one spinor-out, like Mq·δ[a,b] + Z·γ[μ,a,b]·vec[μ]). *)
diracNumeratorSumQ[p_Plus] := Module[{terms = List @@ p, opens},
  opens = openSpinorOf /@ terms;
  AllTrue[terms, ! FreeQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac] &] &&
    AllTrue[opens, Length[#] === 2 &] && SameQ @@ (Sort /@ opens)];
(* Collection of DRESSED vertex sums (ntDiracSlot with k>=1 open legs) is OPT-IN via NT_VERTEX_COLLECT:
   the codegen expansion materialises the structure x dressing product (~10^6 per net on full-basis
   ZAAqbq) and can OOM the kernel. It wins only on small flows. *)
$ntVertexCollect := ntEnvFlag["NT_VERTEX_COLLECT"];
(* A Dirac Plus is collected (kept eager as one token) iff it is not a sector bridge and EITHER
   - a dressed propagator numerator (k=0 open legs, ident/slash only) -> ntDressedNum; OR
   - a Dirac slot (shared spinor pair and open-leg set) -> ntDiracSlot, which for a DRESSED sum
     requires $ntVertexCollect, and for an ALL-NUMERIC sum (e.g. a multi-term projector) is always on:
     it has no dressing dimension, so collecting it is bounded, whereas leaving it a raw Plus makes
     splitColourGroups expand its Cartesian product with the rest of the diagram.
   The `=!= $Failed` decompositions are the exact tests; the other predicates are cheap pre-filters. *)
collectibleDiracSumQRaw[p_Plus] := ! sectorBridgeQ[p] &&
  ((dressedStructureSumQ[p] && diracNumeratorSumQ[p] && dressedNumDecompose[p] =!= $Failed) ||
   (TrueQ[$ntVertexCollect] || ! dressedStructureSumQ[p]) &&
     diracSlotSumQ[p] && diracSlotDecompose[p] =!= $Failed);
collectibleDiracSumQ[_] := False;
distributeQRaw[p_] := sectorBridgeQ[p] ||
  (dressedStructureSumQ[p] && ! (TrueQ[$ntDressCollect] && collectibleDiracSumQ[p]));

(* The SET of γ-count parities the diagram's terms would carry IF every EAGER Dirac Plus were
   distributed — {0}, {1}, or {0,1}. A closed trace is identically zero only when the set is exactly
   {1}; a diagram-global γ count would sum ACROSS the branches of an eager Plus (1 + p̸ against two γ's
   counts 3) and report a parity no actual term has. Computed WITHOUT Expand: the sets are capped at two elements and fold pairwise, so this is linear
   in the factor count. ntDressedNum cannot appear yet (rewriteDressedNums runs later). *)
diracParities[e_Plus] := Union @@ (diracParities /@ (List @@ e));
diracParities[e_Times] := Fold[Union[Flatten[Mod[Outer[Plus, #1, #2], 2]]] &,
                               {0}, diracParities /@ (List @@ e)];
(* b^n: n copies of a definite parity p give n*p mod 2; if b has both parities, so does b^n. Do NOT
   recurse on Times @@ ConstantArray[b, n]: Times re-collapses it to b^n (infinite recursion), and
   this fires on ordinary scalar powers too. *)
diracParities[Power[b_, n_Integer?Positive]] := With[{s = diracParities[b]},
  If[Length[s] > 1, {0, 1}, {Mod[n * First[s], 2]}]];
(* {0, Infinity}, NOT Infinity: this catch-all can receive a BARE ntGamma head, which level spec
   Infinity (= {1, Infinity}) would skip.
   Counts _ntGamma ONLY: ntSigma (2 γ), ntDeltaDirac, ntC (= γ^2 γ^4) and γ5 are all parity-even
   (block-diagonal in the Weyl basis). Keep in step with the C++ counter nAntidiag
   (numeric_contract.hpp / network/dirac.hpp). *)
diracParities[e_] := {Mod[Count[e, _ntGamma, {0, Infinity}], 2]};

(* The odd-trace verdict itself: a γ5-free diagram every one of whose branches is an odd closed
   γ trace is identically zero and can be pruned before any trace runs. *)
vanishingOddTraceQ[diagram_] := FreeQ[diagram, _ntGamma5] && diracParities[diagram] === {1};

(* The slash momentum of a γ·p term: the {coeff, q} pairs of the ntVec factors sharing the γ's Lorentz
   index μ (a signed linear combination like p−l, exactly like ntSigmaLeg). $Failed if ill-formed. *)
dressedVlc[facs_, mu_] := Module[{lf = Times @@ Select[facs, ! FreeQ[#, ntVec[_, mu]] &], terms},
  If[lf === 1, Return[$Failed]];
  terms = If[Head[lf] === Plus, List @@ lf, {lf}];
  Catch[Function[trm, Module[{q = Cases[trm, ntVec[qq_, m_] /; m === mu :> qq, {0, Infinity}],
                              c = trm /. ntVec[_, m_] /; m === mu :> 1},
     If[Length[q] =!= 1 || ! FreeQ[c, ntVec] || ! NumericQ[c], Throw[$Failed, "dvlc"]];
     {c, First[q]}]] /@ terms, "dvlc"]];

(* Decompose one term of a dressed numerator into {din, dout, scalar, spec, otherTensors}: the Dirac
   in/out spinor labels; the scalar (dressing × numeric × denominator) coefficient; the Dirac structure
   spec — {"ident"} (spinor-δ) or {"slash", {{c,q}…}} (γ·p̸); and the NON-Dirac tensor factors (a colour
   δ on the quark line, etc.) that must be COMMON across the sum's terms and factor out. $Failed for an
   unsupported structure (σ, γ5, >1 γ). *)
dressedNumTerm[t_] := Module[{facs = If[Head[t] === Times, List @@ t, {t}], dt, scal, mu, vlc, other},
  scal = Times @@ Select[facs, scalarQ];
  dt = Select[facs, MatchQ[#, _ntDeltaDirac | _ntGamma] &];
  Which[
    MatchQ[dt, {ntDeltaDirac[_, _]}],
      other = Select[DeleteCases[facs, dt[[1]]], tensorQ];
      {dt[[1, 1]], dt[[1, 2]], scal, {"ident"}, Sort[other]},
    MatchQ[dt, {ntGamma[_, _, _]}],
      mu = dt[[1, 1]]; vlc = dressedVlc[facs, mu];
      If[vlc === $Failed, $Failed,
        (other = Select[DeleteCases[facs, dt[[1]]], tensorQ[#] && FreeQ[#, ntVec[_, mu]] &];
         {dt[[1, 2]], dt[[1, 3]], scal, {"slash", vlc}, Sort[other]})],
    True, $Failed]];

(* Multiset of multiplicative factors common to every list in `factLists` (min multiplicity). *)
commonFactorMultiset[factLists_] := Module[{cnts = Counts /@ factLists, keys},
  keys = Intersection @@ (Keys /@ cnts);
  Flatten[Function[k, ConstantArray[k, Min[(Lookup[#, k, 0] &) /@ cnts]]] /@ keys]];

dressedNumDecomposeRaw[p_Plus] := Module[
   (* Flatten only this local numerator sum, and only over its tensor structure. At finite T,
      psdash[p] contains gamma.mu vecs[p, mu]; after fixed-component normalization the temporal
      subtraction must become a separate slash option rather than making the complete propagator
      numerator non-collectible. A tensor-free factor (a composite dressing such as
      -(Zq[q] + Zq[k] RF/|q|)) stays whole: expanding it too would turn each of its atoms into its
      own dressing slot and multiply the distinct sub-terms (4.75x on za4_147). *)
   {rows = dressedNumTerm /@ (List @@ Expand[p, $ntTensorPat]), din, dout, others, common, scalFacs, commonScal, opts},
  If[MemberQ[rows, $Failed], Return[$Failed]];
  {din, dout} = rows[[1, {1, 2}]];
  If[! AllTrue[rows, #[[1]] === din && #[[2]] === dout &], Return[$Failed]]; (* all terms din→dout *)
  others = rows[[All, 5]];
  If[! AllTrue[others, Sort[#] === Sort[others[[1]]] &], Return[$Failed]];   (* common colour factors *)
  common = others[[1]];
  (* factor the scalar coefficient common to every term (the propagator denominator, a flavour δ, …)
     OUT of the sum so it multiplies the whole numerator; only the per-structure residual (the
     dressing that differs: Mq vs Zq·…) stays inside the ntDressedNum options. *)
  scalFacs = Function[s, If[Head[s] === Times, List @@ s, {s}]][#[[3]]] & /@ rows;
  commonScal = commonFactorMultiset[scalFacs];
  opts = MapThread[Function[{r, sf},
     {Times @@ Fold[DeleteCases[#1, #2, {1}, 1] &, sf, commonScal], r[[4]]}], {rows, scalFacs}];
  (Times @@ common) * (Times @@ commonScal) * ntDressedNum[opts, din, dout]];

(* Rewrite each collectible dressed Dirac numerator (a surviving Plus factor) into one ntDressedNum
   token, and SPLICE the factored-out common colour / flavour / denominator product back into the
   factor list as separate factors — so the scalars (flavour δ, 1/denom) land in the diagram coeff
   (where contractFlavour can collapse a flavour chain) and the colour δ folds as its own factor. A
   factor that decomposes to $Failed is left as-is. Only active under $ntDressCollect. *)
rewriteDressedNums[factors_List] := If[! TrueQ[$ntDressCollect], factors,
  Flatten[Function[f, If[Head[f] === Plus && collectibleDiracSumQ[f],
     (* k=0 propagator numerator → ntDressedNum; else k>=1 vertex → ntDiracSlot *)
     With[{r = With[{r0 = dressedNumDecompose[f]}, If[r0 =!= $Failed, r0, diracSlotDecompose[f]]]},
       If[r === $Failed, {f}, If[Head[r] === Times, List @@ r, {r}]]], {f}]] /@ factors]];

(* ---- general collected Dirac slot (ANY open-leg count) --------------------------------
   A collected Dirac slot is a coefficient-weighted sum of Dirac structures that all share the spinor
   in/out pair AND the SAME SET of open Lorentz legs `{μ...}` (k>=0). k=0 is the propagator numerator
   (handled by the ntDressedNum path); k>=1 is a vertex with k open gluon legs (Aqbq: 1; AAqbq: 2; …).
   The engine closes any number of open legs, so the front end only keeps the sum EAGER and packages
   each structure as one option. Unlike the propagator
   collection, an option's structure is kept WHOLE (its Dirac chain × its Lorentz-net factors, e.g. the
   gluon propagator on the open leg); the codegen backend splits it into Dirac tokens vs net factors. *)

(* colour (SU(N)) labels anywhere under an expression — the axes that must factor out of the slot. *)
colourLabelsOf[e_] := DeleteDuplicates @ Flatten @ Cases[e,
  h_ /; (fundamentalSUNQ[h] || adjointSUNQ[h]) :> labelsOf[h], {0, Infinity}];
(* the OPEN Lorentz legs of a term: free indices that are neither spinor nor colour — the gluon axes. *)
openLorentzOf[t_] := Complement[freeIdx[t], allSpinorLabels[t], colourLabelsOf[t]];

(* A collectible general Dirac slot: a Plus whose EXPANDED terms are each a Dirac structure with the
   same 2 open spinor indices and the same NON-EMPTY set of open Lorentz legs (so the surrounding net
   contracts a fixed leg set for every structure choice). Expand first so a term carrying an inner Dirac
   Plus (e.g. a σ commutator written out) splits into monomials. *)
diracSlotSumQ[p_Plus] := Module[{terms = List @@ Expand[p, $ntTensorPat], opens, lors},
  opens = openSpinorOf /@ terms;
  lors  = Sort /@ (openLorentzOf /@ terms);
  AllTrue[terms, ! FreeQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac] &] &&
    AllTrue[opens, Length[#] === 2 &] && SameQ @@ (Sort /@ opens) &&
    Length[First[lors]] >= 1 && SameQ @@ lors];
diracSlotSumQ[_] := False;

(* Decompose a collectible vertex sum into `(commonColour) (commonScalar) ntDiracSlot[opts, din, dout,
   legs]`, where each option is `{residualScalarCoeff, structureProduct}` and `structureProduct` is the
   term's Dirac + Lorentz-net factors (colour and the common scalar factored out). $Failed if the colour
   factor is not common across terms (then the sum is left to distribute). *)
NumTrace::slotorient = "diracSlotDecompose: a collected Dirac slot has `1` candidate in-legs, not 1. A slot is an open chain din->dout, so exactly one open spinor label must be an IN leg and no OUT leg; two (an anomalous qq vertex) or zero (its qbar qbar conjugate) leave it without an orientation, and guessing one could emit the segment backwards. Open spinor labels: `2`. First term: `3`.";

diracSlotDecomposeRaw[p_Plus] := Module[
  {terms = List @@ Expand[p, $ntTensorPat], legs, opens, din, dout, io, ins, outs, dinCands, rows, cols, common, scals, commonScal, opts},
  If[! diracSlotSumQ[p], Return[$Failed]];
  legs   = Sort @ openLorentzOf[First[terms]];
  opens  = openSpinorOf[First[terms]];
  (* ORIENTED din/dout: each Dirac head is an (in,out) spinor edge; the chain runs din→dout. *)
  io = Cases[If[Head[First[terms]] === Times, List @@ First[terms], {First[terms]}],
        ntGamma[_, a_, b_] | ntGamma5[a_, b_] | ntC[a_, b_] | ntDeltaDirac[a_, b_] |
        ntSigma[_, _, a_, b_] :> {a, b}];
  ins = io[[All, 1]]; outs = io[[All, 2]];
  (* ORIENTATION GUARD: din must be the UNIQUE open label that is some head's `in` and no head's `out`.
     Zero or two candidates (q̄q̄ / qq vertex) have no orientation, and since the tokens are spliced in
     CHAIN ORDER a guess would silently emit the segment backwards. *)
  dinCands = Select[opens, MemberQ[ins, #] && ! MemberQ[outs, #] &];
  If[Length[dinCands] =!= 1,
    Message[NumTrace::slotorient, Length[dinCands], Short[opens, 4], Short[First[terms], 6]]; Abort[]];
  din  = First[dinCands];
  dout = First[DeleteCases[opens, din], Last[opens]];
  (* per term -> {scalar, sorted colour factors, sorted structure (Dirac + Lorentz-net, no colour/scalar)} *)
  rows = Function[t, Module[{facs = If[Head[t] === Times, List @@ t, {t}], scal, col, struct},
     scal   = Times @@ Select[facs, scalarQ];
     col    = Sort @ Select[facs, (fundamentalSUNQ[#] || adjointSUNQ[#]) &];
     struct = Sort @ Select[facs, (tensorQ[#] && ! (fundamentalSUNQ[#] || adjointSUNQ[#])) &];
     {scal, col, struct}]] /@ terms;
  cols = rows[[All, 2]];
  If[! AllTrue[cols, # === cols[[1]] &], Return[$Failed]];   (* colour must factor out of the slot *)
  common = cols[[1]];
  (* factor the scalar common to every term (propagator denominator, flavour δ, …) out; the per-term
     dressing residual (Zqbq1 vs Zqbq4 vs …) stays inside the option. *)
  scals = Function[s, If[Head[s] === Times, List @@ s, {s}]][#[[1]]] & /@ rows;
  commonScal = commonFactorMultiset[scals];
  opts = MapThread[Function[{r, sf},
     {Times @@ Fold[DeleteCases[#1, #2, {1}, 1] &, sf, commonScal], Times @@ r[[3]]}], {rows, scals}];
  (Times @@ common) * (Times @@ commonScal) * ntDiracSlot[opts, din, dout, legs]];
diracSlotDecompose[_] := $Failed;

(* ---- per-call memo for the Plus classifiers ------------------------------------------------------
   FunKit reuses index names across diagrams, so one vertex sum is the SAME expression in every diagram
   that contains it (ZA4_147: 2631 Plus factors, 40 distinct). The four functions below (and
   labelCensus on a Plus) are pure in {p, $ntDressCollect, $ntVertexCollect}, so NumTrace and
   FromFunKit Block $ntPlusMemo to <||> and share the answers for one call; outside such a Block
   nothing is cached. The store happens only after f[p] returns, so an Abort is never memoised. *)
$ntPlusMemo = None;
ntPlusMemo[tag_, f_, p_] :=
  If[$ntPlusMemo === None,
    f[p],
    With[{key = {tag, p, TrueQ[$ntDressCollect], TrueQ[$ntVertexCollect]}},
      Lookup[$ntPlusMemo, Key[key], $ntPlusMemo[key] = f[p]]]];
distributeQ[p_]              := ntPlusMemo["dq", distributeQRaw, p];
collectibleDiracSumQ[p_Plus] := ntPlusMemo["cds", collectibleDiracSumQRaw, p];
dressedNumDecompose[p_Plus]  := ntPlusMemo["dnd", dressedNumDecomposeRaw, p];
diracSlotDecompose[p_Plus]   := ntPlusMemo["sd", diracSlotDecomposeRaw, p];

(* Distribute every colour<->Lorentz-bridging sum into its surrounding product (only that
   sum; single-sector sums are left intact as an eager add(...)). Turns a bridging diagram into a
   small linear sum of sector-separable diagrams. Explicit recursion (not a //. rule over
   Orderless Times, which backtracks catastrophically on a large net). *)
expandBridges[e_Plus] := Plus @@ (expandBridges /@ (List @@ e));
expandBridges[e_Times] := Module[{factors = List @@ e, bridge},
  bridge = FirstCase[factors, p_Plus /; distributeQ[p], Missing[]];
  If[MissingQ[bridge],
    Times @@ (expandBridges /@ factors),
    expandBridges[Plus @@ (Times @@ Append[DeleteCases[factors, bridge, {1}, 1], #] & /@ (List @@ bridge))]]];
(* Times auto-collects two byte-identical bridging sums into Power[sum, 2], which the p_Plus selector
   above never picks, so it would pass through UN-DISTRIBUTED (or be read by compileLorentz as a
   closed self-contraction). Two vertices sharing every label are malformed anyway, and Times would
   re-collapse a naive expansion. No known flow produces this; fail loudly. *)
NumTrace::bridgepow = "expandBridges: a colour<->Lorentz-bridging sum appears raised to the power \
`1`, i.e. as `1` byte-identical factors sharing every index label. Such a sum cannot be distributed \
(and identical labels on two distinct vertices are themselves malformed). Offending base:\n`2`";
expandBridges[e : Power[b_Plus, n_Integer]] /; n >= 2 && distributeQ[b] :=
  (Message[NumTrace::bridgepow, n, Short[b, 6]]; Abort[]);
expandBridges[e_] := e;

(* ---- fixed Lorentz components (the finite-T γ0/γi split) ---------------------

   The four-quark Fierz bases (and any finite-T 3+1 split) pin Lorentz indices to a CONCRETE
   component: ntGamma[0, d1, d2] is γ^0. The label machinery has no such notion (it would read the 0
   as a contraction label), so a fixed component is REWRITTEN as a contraction with the constant unit
   basis vector e_i:

       γ^i          ->  ntGamma[μ, d1, d2] ntVec[ntUnitVec[i], μ]     (μ a fresh private dummy)
       g^{i ν}      ->  ntVec[ntUnitVec[i], ν]
       g^{i j}      ->  δ_ij                                          (Euclidean metric)
       σ's free leg ->  the existing SLASH leg against ntUnitVec[i]

   ntUnitVec[i] is an ordinary MOMENTUM symbol whose frame components NumTrace injects as
   UnitVector[4, i+1], so a fixed-component γ is emitted as an ordinary slash, with no C++ change.

   Applied BEFORE expandBridges/checkLabels, so no integer Lorentz slot ever reaches the label
   machinery and tensorQ/labelsOf/freeIdx/labelCensus stay untouched. (ntVec[q, i_Integer] is NOT
   rewritten: it is already the scalar component q_i, resolved by the frame in the coefficient.) *)

(* Each rewritten γ needs its OWN dummy — Unique[] inside the RHS fires once per match. *)
expandFixedComponents[e_] := e //. {
  ntGamma[i_Integer, d1_, d2_] :>
    With[{mu = Unique["fixc$"]}, ntGamma[mu, d1, d2] ntVec[ntUnitVec[i], mu]],
  ntMetric[i_Integer, j_Integer] :> If[i === j, 1, 0],
  ntMetric[i_Integer, nu_] :> ntVec[ntUnitVec[i], nu],
  ntMetric[nu_, i_Integer] :> ntVec[ntUnitVec[i], nu],
  ntSigma[{"free", i_Integer}, legB_, d1_, d2_] :>
    ntSigma[{"slash", {{1, ntUnitVec[i]}}}, legB, d1, d2],
  ntSigma[legA_, {"free", i_Integer}, d1_, d2_] :>
    ntSigma[legA, {"slash", {{1, ntUnitVec[i]}}}, d1, d2],
  ntEpsilon[a___, i_Integer, b___] :>
    With[{mu = Unique["fixc$"]}, ntEpsilon[a, mu, b] ntVec[ntUnitVec[i], mu]]};

(* The unit vectors a rewritten network needs, as frame entries. Component i (0-based, 0 = the
   temporal/Matsubara direction) is UnitVector[4, i+1]. *)
unitVecFrame[net_] := Association[
  (# -> UnitVector[4, First[#] + 1]) & /@ DeleteDuplicates[Cases[net, _ntUnitVec, {0, Infinity}]]];

(* A component index outside 0..3 is a caller error (a mis-set 3+1 convention), not a tracer bug. *)
NumTrace::fixcomp = "Fixed Lorentz component `1` is out of range: a component index must be 0..3 \
(0 = temporal). Check the basis's 3+1 convention.";

(* ---- finite-T SPATIAL vectors (FormTracer's `vecs`) -------------------------------------------

   vecs[q, mu] is the spatial part of q as a 4-vector, {0, q_1, q_2, q_3}. FromFunKit rewrites it to
   ntVec[ntSpatialVec[q], mu], a new momentum LEAF, so a spatial slash is an ordinary slash (as with
   ntUnitVec above); the zero temporal component is a structural zero in the frame spec. ntSPS, the
   spatial scalar product, needs no leaf; gen_spatialvec_numeric.wls pins that the two agree.

   The projection is linear, so push it through sums and numeric factors FIRST: then only BASE
   momenta become frame keys (no duplicate leaf for ntSpatialVec[p - l]), and the leaf is an atom by
   the time canonicalizeMomentumSigns' negMomQ inspects it. *)
expandSpatialVecs[e_] := e //. {
  ntSpatialVec[ntSpatialVec[q_]]      :> ntSpatialVec[q],   (* idempotent: the bar of a bar *)
  ntSpatialVec[0]                     :> 0,
  ntSpatialVec[a_Plus]                :> (ntSpatialVec /@ a),
  ntSpatialVec[c_?NumericQ * q_]      :> c ntSpatialVec[q],
  ntSpatialVec[c_?NumericQ]           :> 0};

(* An ntSpatialVec whose argument the frame cannot resolve would silently become a leaf with
   SYMBOLIC components — the momentum symbol itself sitting in a component slot — which survives all
   the way into the emitted arithmetic. Refuse it here instead. *)
NumTrace::spatialframe = "ntSpatialVec[`1`]: the frame does not resolve `1` to four components \
(got `2`). Every momentum appearing under a spatial vector (FormTracer's vecs[q, mu]) must be a \
frame key or a linear combination of frame keys.";

(* The spatial vectors a rewritten network needs, as frame entries: the parent momentum's components
   with the temporal slot zeroed. Call AFTER unitVecFrame has joined the frame, so that a spatial
   vector of a unit basis vector resolves too.

   These entries make the spatial vector a frame/env citizen (buildEnv gives it a Base). The numeric
   backend does NOT compute with them (mutating them is byte-identical-inert): its component table
   comes from the frame SPEC (CodegenFrames.m unitLoopFrameSpec / unitLoopMixedFrameSpec /
   polyFrameSpec), which derives the spatial vector from its parent the same way. Keep the two
   derivations in step. *)
spatialVecFrame[net_, frame_] := Association[
  Function[sv, Module[{c = resolveComponents[First[sv], frame]},
      If[! MatchQ[c, {_, _, _, _}],
        Message[NumTrace::spatialframe, First[sv], c]; Abort[]];
      sv -> ReplacePart[c, 1 -> 0]]] /@
    DeleteDuplicates[Cases[net, _ntSpatialVec, {0, Infinity}]]];

(* ---- SU(N) FUNDAMENTAL Levi-Civita -------------------------------------------

   ntEpsFund[N, i1, ..., iN] (exactly N indices) comes from FunKit's epsFundCol/epsFundFlav (the
   four-quark Fierz bases' diquark channels). It is a REWRITE, not an engine token: an epsilon is
   contracted ONLY in pairs (a lone one is not an SU(N) invariant), and a pair folds to a determinant
   of Kronecker deltas,

       eps_{a1..ak c1..cm} eps_{a1..ak d1..dm}  =  k! * det( delta_{c_p d_q} )        (m = N - k)

   i.e. a Plus of ntSUNDeltaFund products with numeric coefficients, which lands on the
   constant-colour branch-list path (compileColourSum in CodegenNets.m).
   Applied BEFORE expandBridges/checkLabels so the object those validate is the one that compiles.
   ntEpsFund is "Rewritten" in $ntHeads and so NOT in the codegen colour tables ($colourHeadPat /
   colourFacStr / labelDimAssoc): a survivor of this rewrite fails (MakeNTKernel::colleak) instead of emitting. *)

(* The dimension of the index space an epsilon head lives in. *)
epsDimOf[ntEpsFund[n_, __]] := n;

$ntEpsMaxPairTerms = 720;   (* 6!; 7! = 5040 would breach $ntColSumMaxBranches (4096) on its own *)

NumTrace::spinorbase = "A diagram carries `1` Lorentz/colour axis labels, which reaches the spinor axis id base `2`. Axis ids are how the engine decides what contracts with what, so the two ranges meeting means a Lorentz axis and a spinor axis share an id and get fused into one contraction — a silently WRONG number, not an error. Raise NumTracer`Private`$ntSpinorIdBase above the label count (the ids are dictionary keys, so a larger base costs nothing). Aborting instead.";

(* First axis id handed to a SPINOR (Dirac) label. Lorentz/colour labels are numbered from 0 upward,
   so this caps their count per diagram (analyseDiagram aborts before the ranges meet). Real diagrams
   peak in the low tens; the ids are dictionary keys, so a large base is free. *)
$ntSpinorIdBase = 100;
NumTrace::epsbig = "expandFundEps: an epsilon pair in dimension `1` sharing `2` index/indices \
expands to `3`! = `4` Kronecker-delta terms (limit `5`). Emitting these would hand the colour \
branch-list lowering a list it would either reject far downstream with an opaque branch count, or \
— worse — accept and turn into a generator source large enough to OOM the C++ compiler, with \
nothing pointing at a Levi-Civita as the cause. Contract more indices between the two epsilons, or \
raise $ntEpsMaxPairTerms if a flow genuinely needs this.";

(* The pair contraction, written DIMENSION-PARAMETRIC (dim and the delta constructor), so it is
   unit-testable at dim = 2..5 against LeviCivitaTensor. *)
epsPairExpand[dim_Integer, uu_List, vv_List, deltaOf_] := Module[
  {shared, cA, cB, posA, posB, sgn, k, m},
  shared = Intersection[uu, vv];
  k = Length[shared]; m = dim - k;
  If[m! > $ntEpsMaxPairTerms,
    Message[NumTrace::epsbig, dim, k, m, m!, $ntEpsMaxPairTerms]; Abort[]];
  cA = DeleteCases[uu, Alternatives @@ shared];   (* order-preserving complements *)
  cB = DeleteCases[vv, Alternatives @@ shared];
  (* the sign of moving the shared labels to the front of each epsilon, computed from POSITIONS (a
     Signature on the symbol list itself would sort by symbol NAME, which is meaningless here) *)
  posA = Flatten[Position[uu, #, {1}, 1] & /@ Join[shared, cA]];
  posB = Flatten[Position[vv, #, {1}, 1] & /@ Join[shared, cB]];
  sgn = Signature[posA] Signature[posB];
  sgn k! Total[(Signature[#] Times @@ MapThread[deltaOf, {cA, cB[[#]]}]) & /@ Permutations[Range[m]]]];

NumTrace::epsrank = "expandFundEps: a fundamental Levi-Civita with `1` indices at rank N = `2`. The \
fundamental epsilon of SU(N) carries EXACTLY N indices (SU(3) colour: 3; SU(2) isospin: 2). A \
mismatch means the rank injected by FromFunKit disagrees with the basis (e.g. a flavour epsilon \
routed to the colour rank, or Nc unset). Contracting it anyway would build a determinant of the \
wrong size and return a silently wrong number. Offending factor:\n`3`";
NumTrace::epsodd = "expandFundEps: `1` FUNDAMENTAL Levi-Civita factor(s) survived the pair \
contraction. A fundamental epsilon is contracted ONLY in pairs (eps.eps = k! det(delta)); a lone \
one is not an SU(N) invariant and has no representation in the engine. The input carries an odd \
number of them AT ONE RANK — a basis/contraction error. (A pair STRADDLING an eager Plus is not \
this case: it is joined by distributing into the sum, see expandFundEpsRec. The ADJOINT epsilon \
does not come through here at all — at rank 2 it is rewritten to f^abc in FromFunKit.) \
Offending factor(s):\n`2`";
NumTrace::epsambig = "expandFundEps: `1` fundamental Levi-Civita factors of the SAME rank N = `2` \
remain, and the best available pairing shares NO index, so which two are partners is not \
determined. Rank alone does not identify the index SPACE — at Nc == Nf a colour and a flavour \
epsilon are indistinguishable here, and pairing across the two spaces is not an identity at all: \
it contracts colour indices against flavour ones and returns a silently wrong number. Refusing \
rather than guessing. Offending factors:\n`3`";

(* Contract every epsilon pair, one multiplicative context at a time. *)
expandFundEps[e_] /; FreeQ[e, _ntEpsFund] := e;
expandFundEps[e_] := Module[{res},
  res = expandFundEpsRec[e];
  With[{left = Cases[res, _ntEpsFund, {0, Infinity}]},
    If[left =!= {}, Message[NumTrace::epsodd, Length[left], Short[left, 4]]; Abort[]]];
  res];

expandFundEpsRec[e_Plus] := expandFundEpsRec /@ e;
(* eps^2 is a self-contraction of an epsilon with itself: every index is shared, so k = N, m = 0 and
   the value is N!. Compute it ARITHMETICALLY: expandFundEpsRec[b b] would re-collapse to b^2 and
   recurse forever. Higher powers fall through to the epsodd guard. *)
expandFundEpsRec[Power[b_ntEpsFund, 2]] := Module[{idx = Rest[List @@ b], d = epsDimOf[b]},
  If[Length[idx] =!= d, Message[NumTrace::epsrank, Length[idx], d, b]; Abort[]];
  If[Length[DeleteDuplicates[idx]] =!= Length[idx], 0,
     epsPairExpand[d, idx, idx, ntSUNDeltaFund[d, #1, #2] &]]];
expandFundEpsRec[e_Times] := Module[{fs, eps, rest, cand, pair, uu, vv, plusEps, host, others},
  (* recurse into the factors FIRST: an epsilon pair frequently lives inside a Plus factor (a
     multi-term projector). *)
  fs = expandFundEpsRec /@ (List @@ e);
  eps = Cases[fs, _ntEpsFund];
  (* Times @@ fs, NOT e, so the recursion's work is kept. FreeQ over all of fs: an epsilon inside a
     Plus factor may still await a partner from out here (the straddle case at the bottom). *)
  If[FreeQ[fs, _ntEpsFund], Return[Times @@ fs, Module]];
  (* validate arity, and kill a degenerate epsilon (a repeated index) before anything else *)
  Function[h, With[{idx = Rest[List @@ h]},
     If[Length[idx] =!= epsDimOf[h],
       Message[NumTrace::epsrank, Length[idx], epsDimOf[h], h]; Abort[]];
     If[Length[DeleteDuplicates[idx]] =!= Length[idx], Return[0, Module]]]] /@ eps;
  rest = DeleteCases[fs, _ntEpsFund];
  (* Pair greedily by shared-index count, highest first: the pair sharing most indices contracts to
     the fewest terms ((N-k)!). Any pairing WITHIN ONE INDEX SPACE gives the same VALUE — eps.eps =
     k! det(delta) is an identity that holds whatever else multiplies it — so only the term count
     depends on the choice.
     Candidates are restricted to EQUAL RANK: across ranks epsPairExpand would silently drop the
     surplus indices of one partner. Unequal ranks are left unpaired for the epsodd guard. *)
  While[Length[eps] >= 2,
    cand = Select[Subsets[Range[Length[eps]], {2}],
             epsDimOf[eps[[#[[1]]]]] === epsDimOf[eps[[#[[2]]]]] &];
    If[cand === {}, Break[]];
    pair = First@MaximalBy[cand,
             Length[Intersection[Rest[List @@ eps[[#[[1]]]]], Rest[List @@ eps[[#[[2]]]]]]] &];
    uu = Rest[List @@ eps[[pair[[1]]]]]; vv = Rest[List @@ eps[[pair[[2]]]]];
    (* Equal rank is necessary but NOT sufficient to identify partners: at Nc == Nf a colour and a
       flavour epsilon carry the same rank. Sharing an index proves they meet; sharing none leaves it
       undetermined, so a zero-overlap choice is only safe when it is FORCED (exactly two of this
       rank left). Otherwise refuse. *)
    If[Intersection[uu, vv] === {} &&
       Count[eps, h_ /; epsDimOf[h] === epsDimOf[eps[[pair[[1]]]]]] > 2,
      Message[NumTrace::epsambig,
        Count[eps, h_ /; epsDimOf[h] === epsDimOf[eps[[pair[[1]]]]]],
        epsDimOf[eps[[pair[[1]]]]], Short[Select[eps, epsDimOf[#] === epsDimOf[eps[[pair[[1]]]]] &], 6]];
      Abort[]];
    AppendTo[rest, epsPairExpand[epsDimOf[eps[[pair[[1]]]]], uu, vv,
                     With[{n = epsDimOf[eps[[pair[[1]]]]]}, ntSUNDeltaFund[n, #1, #2] &]]];
    eps = Delete[eps, {{pair[[1]]}, {pair[[2]]}}]];
  (* STRADDLE: an epsilon's partner may sit inside an eager Plus factor (or in a DIFFERENT Plus
     factor) rather than out here, so no amount of pairing at one multiplicative level can join
     them — e.g. eps_col[a,A1,A3] eps_flav[F1,F3] * (eps_col[a,A2,A4] eps_flav[F2,F4] D1 - ...),
     the four-quark diquark vertex. Distribute the remaining factors into ONE such sum and recurse;
     each step removes one epsilon-bearing Plus, so this terminates. Times, never Expand, so the
     summands' own internal sums stay intact. *)
  If[eps =!= {} || ! FreeQ[rest, _ntEpsFund],
    plusEps = Select[rest, Head[#] === Plus && ! FreeQ[#, _ntEpsFund] &];
    If[plusEps =!= {},
      host   = First[plusEps];
      others = Times @@ Join[DeleteCases[rest, host, {1}, 1], eps];
      Return[Plus @@ (expandFundEpsRec[others #] & /@ (List @@ host)), Module]]];
  Times @@ Join[rest, eps]];   (* a leftover odd epsilon rides along to the epsodd guard *)
expandFundEpsRec[e_] := e;

(* ---- momentum sign canonicalisation ----------------------------------------- *)

(* buildEnv keys momentum Bases on the momentum EXPRESSION, so `q` and `-q` would get separate Bases
   and Inv slots, and the sign twins would block the sub-term and emitted-body dedup (both key on the
   base name). Every momentum-carrying head is LINEAR in its momentum (vec(-q) = -vec(q)) or EVEN
   (the projectors, 1/q²), so a canonical sign is exact: flip iff the coefficient of the first
   variable (Sort order) is negative.
   NT_NO_SIGN_CANON: escape hatch back to expression-keyed bases (used by the lambda3d_small control). *)

negMomQ[q_] := Module[{vars = Sort[Variables[q]], c},
  vars =!= {} && (c = Coefficient[q, First[vars]]; NumericQ[c] && c < 0)];

canonSlashPairs[vlc_List] :=
  Replace[vlc, {c_, q_} /; negMomQ[q] :> {-c, Expand[-q]}, {1}];

(* One rule per momentum-carrying head: q -> -q flips the head by its registered sign. *)
$ntMomentumSignRules = Function[rec,
    With[{h = ntHeadSym[rec], s = rec["Momentum"],
          k = Replace[Extract[rec, Key["Form"], Hold], Hold[_[args___]] :> Length[Hold[args]]] - 1},
      h[q_, rest : Repeated[_, {k}]] /; negMomQ[q] :> s h[Expand[-q], rest]]] /@
  Select[$ntHeads, KeyExistsQ[#, "Momentum"] &];

canonicalizeMomentumSigns[net_] :=
  If[ntEnvFlag["NT_NO_SIGN_CANON"],
    net,
    net /. {
      Sequence @@ $ntMomentumSignRules,
      (* sigma slash legs carry {coeff, q} pairs directly (not ntVec factors) — linear likewise *)
      ntSigma[{"slash", vlcA_List}, legB_, d1_, d2_] /; AnyTrue[vlcA, negMomQ[Last[#]]&] :>
        ntSigma[{"slash", canonSlashPairs[vlcA]}, legB, d1, d2],
      ntSigma[legA_, {"slash", vlcB_List}, d1_, d2_] /; AnyTrue[vlcB, negMomQ[Last[#]]&] :>
        ntSigma[legA, {"slash", canonSlashPairs[vlcB]}, d1, d2]
    }];

(* ---- env-id layout ---------------------------------------------------------- *)

(* Assign each distinct momentum a Base (4 consecutive Var ids) and, where a
   projector needs it, an Inv id holding 1/q^2. Var ids index the runtime renv[]. *)
buildEnv[momenta_List, invMomenta_List, invSMomenta_List] := Module[{env = <||>, base = 0, inv},
  Do[env[q] = <|"Base" -> base, "Inv" -> None, "InvS" -> None|>; base += 4, {q, momenta}];
  inv = base;
  Do[env[q]["Inv"] = inv++, {q, invMomenta}];      (* full 1/q² slots *)
  Do[env[q]["InvS"] = inv++, {q, invSMomenta}];    (* spatial 1/|q⃗|² slots (finite-T E/M projectors) *)
  {env, inv}  (* inv = total env size NEnv *)
];

(* Per-momentum component mask: bit i set <=> component i is structurally nonzero in the frame. *)
frameMask[components_List] := FromDigits[Reverse[Boole[# =!= 0 && # =!= 0.] & /@ components], 2];

(* ---- NumTrace --------------------------------------------------------------- *)

(* NumTrace is deliberately serial: mapping analyseDiagram / labelCensus over Wolfram subkernels is
   marshalling-bound (55 s serial vs 55 s parallel on full-basis ZAAqbq). *)

(* Each SU(N) group head carries its own rank N as the first argument (baked in when the
   network is built — Global`Nc for colour, the FromFunKit "FlavourGroup" option for the
   isospin group), so NumTrace itself takes no group option. *)
Options[NumTrace] = {"Frame" -> <||>, "Args" -> {}, "Dressings" -> {}, "DressingCollection" -> True};

(* the "DressingCollection" value of the most recent FromFunKit call; None = no FromFunKit yet *)
$ntFromFunKitDressCollect = None;
NumTrace::sunrank = "every SU(N) head must carry an integer rank N >= 1 as its first argument (got `1`). Set Nc (SetNc[3]) before tracing colour heads, and pass \"FlavourGroup\" -> n to FromFunKit for the isospin group.";
NumTrace::nfsym = "the flavour count Nf is not a defined integer (a symbolic Nf is present in the network). Call SetNf[2] (TensorBases) before generating.";
NumTrace::dresscollect = "NumTrace runs with \"DressingCollection\" -> `1`, but the most recent FromFunKit \
call used `2`. FromFunKit already kept (or distributed) the dressed Dirac numerators according to its \
own setting; pass the same value to both.";

NumTrace[net_, OptionsPattern[]] := Block[{$ntProfOn = TrueQ[$NumTracerVerbose], $ntProf = <||>, $ntPlusMemo = <||>}, Module[
  {frame, args, dress, badRanks, net2, diagrams, allMom, invMom, invSMom, env, nenv, diags, ntT0},
  (* whole-NumTrace wall clock, so a gap between the total and the [prof] parts is visible *)
  ntT0   = AbsoluteTime[];
  frame  = OptionValue["Frame"];
  args   = OptionValue["Args"];
  dress  = OptionValue["Dressings"];
  (* symbolic dressing collection (see $ntDressCollect); set here so expandBridges and analyseDiagram
     both see it *)
  $ntDressCollect = TrueQ[OptionValue["DressingCollection"]];
  If[BooleanQ[$ntFromFunKitDressCollect] && $ntFromFunKitDressCollect =!= $ntDressCollect,
    Message[NumTrace::dresscollect, $ntDressCollect, $ntFromFunKitDressCollect]];

  (* SU(N) ranks must be compile-time integers (they pick the correct-dimension typed-out
     group matrices). Every group head's leading argument N is checked up front, so an
     undefined symbol aborts with a clear message rather than being baked into the kernel. *)
  badRanks = DeleteDuplicates @ Select[
    Cases[net, h : $ntSUNHeadPat :> sunRankOf[h], Infinity],
    ! (IntegerQ[#] && # >= 1) &];
  If[badRanks =!= {},
    Message[NumTrace::sunrank, badRanks]; Abort[]];
  (* A closed quark flavour loop folds to Global`Nf; it is undefined only if Nf is still a symbol
     AND appears in the network. *)
  If[! IntegerQ[Global`Nf] && ! FreeQ[net, Global`Nf],
    Message[NumTrace::nfsym]; Abort[]];

  (* Fixed Lorentz components (γ^0 & co) are rewritten FIRST, so no integer Lorentz slot reaches the
     label machinery (see expandFixedComponents); the unit vectors join the frame as momenta. *)
  net2 = ntProfTimed["canonicalizeMomentumSigns", canonicalizeMomentumSigns @
           ntProfTimed["expandFundEps", expandFundEps @
             ntProfTimed["expandSpatialVecs", expandSpatialVecs @
               ntProfTimed["expandFixedComponents", expandFixedComponents[net]]]]];
  With[{bad = DeleteDuplicates @ Cases[net2, ntUnitVec[i_] :> i, {0, Infinity}]},
    If[! AllTrue[bad, IntegerQ[#] && 0 <= # <= 3 &],
      Message[NumTrace::fixcomp, Select[bad, ! (IntegerQ[#] && 0 <= # <= 3) &]]; Abort[]]];
  frame = Join[frame, unitVecFrame[net2]];
  (* FINITE-T SPATIAL VECTORS (FormTracer's vecs) join the frame the same way, but AFTER the unit
     vectors — a spatial vector's components are read off its parent's, so the parent must already
     resolve. See expandSpatialVecs / spatialVecFrame. *)
  frame = Join[frame, spatialVecFrame[net2, frame]];

  (* the top-level sum is the (linear) sum of DIAGRAMS; single-sector vertex sums stay eager, while
     colour<->Lorentz-bridging sums are distributed (see sectorBridgeQ). *)
  With[{ntT = First@AbsoluteTiming[
  diagrams = With[{ex = expandBridges[net2]}, If[Head[ex] === Plus, List @@ ex, {ex}]];]},
    ntLog["[prof] NumTrace expandBridges: ", ntT, " s"]];
  (* odd-trace vanishing: a γ5-free closed Dirac trace of an ODD number of gammas is zero, so drop
     such diagrams. The verdict is PER BRANCH (diracParities): an eager Dirac Plus (a collected
     propagator Mq·δ + Z·p̸, a multi-term projector) can mix parities, and a diagram-global count
     would drop non-vanishing even branches. Only all-odd diagrams are dropped; this is pure pruning,
     since the per-branch trace zeroes odd branches on its own. *)
  diagrams = ntProfTimed["oddTracePrune", Select[diagrams, ! vanishingOddTraceQ[#] &]];

  (* Validate every distributed diagram BEFORE analyseDiagram assigns axis ids (one id per DISTINCT
     label), per diagram because only after expandBridges is "1 = free, 2 = contracted" the
     invariant. See the label-validation section below. *)
  With[{ntT = First@AbsoluteTiming[
  With[{frees = If[TrueQ[$ntCheckLabels],
      (* the census (labelCensus) is pure; the abort/Message validation is kept separate so a
         failure reports the diagram index (see checkLabels) *)
      With[{census = ntProfTimed["labelCensus", labelCensus /@ diagrams]},
        ntProfTimed["checkLabels", MapThread[checkLabels[#1, #2, #3] &, {diagrams, census, Range[Length[diagrams]]}]]],
      ConstantArray[{}, Length[diagrams]]]},
    ntLog["[labels] ", Length[diagrams], " diagram(s) validated; free-index set(s) = ",
      DeleteDuplicates[Sort /@ frees]]];]},
    ntLog["[prof] NumTrace checkLabels: ", ntT, " s"]];

  (* global env layout: every distinct momentum, and which ones need a 1/q^2 slot. Scan net2, not net:
     the unit and spatial vectors introduced above are momenta that need an env Base. Deduplicating
     the heads first keeps each momentum's first-occurrence order. *)
  ntProfTimed["buildEnv", With[{tens = DeleteDuplicates @ Cases[net2, _?tensorQ, Infinity]},
    allMom = DeleteDuplicates[momentumOf /@ tens] // DeleteCases[None];
    invMom = DeleteDuplicates[momentumOf /@ Select[tens, needsInvQ]];
    invSMom = DeleteDuplicates[momentumOf /@ Select[tens, needsInvSQ]];
    {env, nenv} = buildEnv[allMom, invMom, invSMom]]];

  (* The text of this log line is a CONTRACT: tests/gen/regen_check.sh's flow_counts() parses it. *)
  With[{ntT = First@AbsoluteTiming[diags = analyseDiagram /@ diagrams]},
    ntLog["[prof] NumTrace analyseDiagram (", Length[diagrams], " diagrams): ", ntT, " s"]];

  (* NO FLAVOUR DELTA MAY LEAVE HERE: both contractFlavour passes and promoteFlavResidue have run,
     and a surviving flavDelta would be classified as a scalar and leak into the C++ (see
     NumTrace::flavleak). $ntCppLeakPatterns (CodegenCommon.m) is the textual backstop. *)
  With[{leak = DeleteDuplicates @ Cases[diags, _flavDelta, {0, Infinity}]},
    If[leak =!= {}, Message[NumTrace::flavleak, Short[leak, 8]]; Abort[]]];

  ntProfReport["[prof]   NumTrace part "];
  ntLog["[prof] NumTrace TOTAL (", Length[diags], " diagrams): ", AbsoluteTime[] - ntT0, " s"];

  NTKernel[<|
    "Diagrams"  -> diags,
    "Env"       -> env,
    "NEnv"      -> nenv,
    "Frame"     -> frame,
    "Args"      -> args,
    "Dressings" -> dress
  |>]
]];

(* One diagram -> {pure-scalar coeff, axis-id map, tensor components}. Components keep
   their factors un-expanded (heads, Plus-vertices, Times-structures); the recursive
   net builder (CodegenNets.m compileLorentz) turns Plus -> add(...), Times -> mul(...). *)
analyseDiagram[diagram_] := Module[{factors, tensorF, ids},
  factors = ntProfTimed["rewriteDressedNums", rewriteDressedNums @
    ntProfTimed["splitSelfTraces", splitSelfTraces[If[Head[diagram] === Times, List @@ diagram, {diagram}]]]];
  (* Promote unclosable flavour deltas into the SU(N) engine: after rewriteDressedNums (which flattens
     straddling chains) and before the axis-id partition below. A no-op when every line closes. *)
  factors = ntProfTimed["promoteFlavResidue", promoteFlavResidue[factors]];
  tensorF = Select[factors, ! scalarQ[#] &];
  (* Partition labels by sector: spinor axes get ids >= $ntSpinorIdBase, so the engine (which
     contracts by MATCHING ID) never fuses a spinor axis with a Lorentz/colour one. This holds only
     while the non-spinor count stays below the base, hence the check. *)
  ids     = With[{labs = DeleteDuplicates @ Flatten[allLabels /@ tensorF],
                  spn  = DeleteDuplicates @ Flatten[allSpinorLabels /@ tensorF]},
              With[{nonsp = DeleteCases[labs, Alternatives @@ spn]},
                If[Length[nonsp] > $ntSpinorIdBase,
                  Message[NumTrace::spinorbase, Length[nonsp], $ntSpinorIdBase]; Abort[]];
                Join[AssociationThread[nonsp -> Range[0, Length[nonsp] - 1]],
                     AssociationThread[spn -> Range[$ntSpinorIdBase, $ntSpinorIdBase - 1 + Length[spn]]]]]];
  <|
    (* second contractFlavour pass: closes flavour-δ chains the dressing collection factored out of an
       eager numerator (FromFunKit's pass could not see them) *)
    "Coeff"      -> If[TrueQ[$ntDressCollect], contractFlavour[Times @@ Select[factors, scalarQ]],
                       Times @@ Select[factors, scalarQ]],
    "Ids"        -> ids,
    (* "Constant" means "a constant SU(N) component": such a component is handed to compileColour,
       which understands only group heads. So EVERY non-SU(N) tensor head makes a component
       non-constant, not just the momentum-carrying ones: a momentum-free Dirac δ-loop or closed metric loop would
       otherwise be emitted as raw Mathematica. The Dirac case relies on compileDirac restoring
       tr(1) = 4 for token-free loops (CodegenNets.m). *)
    "Components" -> (<|"Factors" -> ntProfTimed["orderFactors", orderFactors[#]],
                       "Constant" -> FreeQ[#, $ntNonSUNHeadPat]|> &
                     /@ ntProfTimed["connectedComponents", connectedComponents[tensorF]])
  |>
];

(* ---- per-diagram label validation ------------------------------------------------
   Run on ONE diagram AFTER expandBridges (a flat Times whose factors are heads, eager
   single-sector Plus vertices, Powers and scalars). The invariant, per diagram:
     * a label occurring ONCE  is a free (external) index of the whole trace;
     * a label occurring TWICE is contracted — the engine pairs the two axes by id;
     * a label occurring 3+ times is MALFORMED: axes sharing an id are silently mis-paired into
       a wrong number, so this must Abort rather than Message.
   An eager Plus is counted ONCE, via the free set its summands must agree on (the add(...)
   alignment precondition freeIdx[_Plus] assumes). A summand's internal dummies are private and
   must not appear anywhere else. *)

NumTrace::flavleak = "a fundamental-flavour Kronecker delta survived BOTH contractFlavour passes \
and the promotion into the SU(N) engine. It is now neither contracted nor a tensor: scalarQ is a \
FreeQ over the nt* heads, so it counts as a scalar COEFFICIENT, its indices are invisible to \
labelsOf/freeIdx/checkLabels, and CForm would print it into the kernel as \
NumTracer_Private_flavDelta(F1, F2) — an undeclared identifier that GCC and Clang both parse, so \
the failure would surface as a link/compile error far from here, or compile silently wrong. Most \
likely cause: a flavour delta buried inside an eager (un-distributed) Plus that rewriteDressedNums \
did not lift out, so promoteFlavResidue saw it as already closed. Offending delta(s):\n`1`";

NumTrace::badlabel = "Diagram `1`: index label `2` occurs `3` times (expected 1 = free, \
2 = contracted). The engine contracts axes by matching id, so `3` axes sharing this label \
are silently mis-paired into a wrong number. Offending diagram:\n`4`";
NumTrace::plusfree = "Diagram `1`: the summands of an eager (un-distributed) sum expose \
DIFFERENT free indices `2` — the eager add(...) cannot align them. Offending sum:\n`3`";
NumTrace::privclash = "Diagram `1`: label(s) `2` are private dummies of one factor but also \
occur outside it — a dummy-name collision between two independently generated objects. \
Offending diagram:\n`3`";

(* {exposed-multiset, private-set, bad-list} of a (sub)expression. *)
labelCensus[p_Plus] := ntPlusMemo["lc", labelCensusPlus, p];
labelCensus[e_] := Which[
  tensorQ[e],
    With[{ls = labelsOf[e]},
      (* a label repeated WITHIN one head is a legal self-trace; splitSelfTraces resolves it *)
      {DeleteDuplicates[ls], Cases[Tally[ls], {l_, c_} /; c >= 2 :> l], {}}],

  Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 1 && ! scalarQ[e],
    (* n copies sharing the SAME labels: a closed self-contraction (see compileLorentz) *)
    Module[{c = labelCensus[e[[1]]], n = e[[2]]},
      Which[
        n === 1, c,
        n === 2, {{}, Union[c[[1]], c[[2]]], c[[3]]},
        True,    {{}, Union[c[[1]], c[[2]]], Join[c[[3]], ({#, n} &) /@ c[[1]]]}]],

  Head[e] === Times,
    Module[{sub = labelCensus /@ Select[List @@ e, ! scalarQ[#] &], exp, tal, priv, bad},
      If[sub === {}, Return[{{}, {}, {}}]];
      exp  = Join @@ sub[[All, 1]];
      tal  = Tally[exp];
      priv = Union @@ sub[[All, 2]];
      bad  = Join[Join @@ sub[[All, 3]], Cases[tal, {l_, c_} /; c > 2 :> {l, c}]];
      (* a child's private dummy must not be exposed or private anywhere else *)
      Do[With[{mine = sub[[i, 2]],
               others = Union[Join @@ Delete[sub, i][[All, 1]], Join @@ Delete[sub, i][[All, 2]]]},
           With[{clash = Intersection[mine, others]},
             If[clash =!= {}, bad = Join[bad, ({#, "private-clash"} &) /@ clash]]]],
         {i, Length[sub]}];
      {Cases[tal, {l_, c_} /; c == 1 :> l],
       Union[priv, Cases[tal, {l_, c_} /; c == 2 :> l]],
       bad}],

  True, {{}, {}, {}}];

labelCensusPlus[e_Plus] := Module[{sub = labelCensus /@ (List @@ e), frees},
  frees = Sort /@ sub[[All, 1]];
  {If[frees === {}, {}, First[frees]],
   Union @@ sub[[All, 2]],
   Join[Join @@ sub[[All, 3]],
        If[Length[DeleteDuplicates[frees]] > 1, {{frees, "plus-free-mismatch"}}, {}]]}];

(* Escape hatch: NT_NO_LABEL_CHECK=1 disables the label census; anything falsy leaves it ON. *)
$ntCheckLabels := !ntEnvFlag["NT_NO_LABEL_CHECK"];

(* Validate a PRECOMPUTED census and return the diagram's free-index set. Split from labelCensus so
   the pure counting stays free of side effects and this half owns the diagnostics; `diagram` is
   carried only for the Short[...] in the message, and `idx` is the diagram's 1-based position (what
   the user sees in the abort). *)
checkLabels[diagram_, census_, idx_] := (
  With[{bad = census[[3]]},
    If[bad =!= {},
      Do[Switch[b[[2]],
           "plus-free-mismatch", Message[NumTrace::plusfree, idx, b[[1]], Short[diagram, 8]],
           "private-clash",      Message[NumTrace::privclash, idx, b[[1]], Short[diagram, 8]],
           _,                    Message[NumTrace::badlabel, idx, b[[1]], b[[2]], Short[diagram, 8]]],
         {b, bad}];
      Abort[]]];
  census[[1]]);
