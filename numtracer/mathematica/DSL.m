(* ::Package:: *)

(* DSL analysis: classify the heads, split each term into independent contraction
   components, and allocate the env-id layout the `et` builders consume.

   Loaded inside NumTracer`Private` by NumTracer.m — public symbols (NumTrace, the
   nt* heads) already exist in the NumTracer` context. *)

(* ---- environment flags ------------------------------------------------------- *)

(* The ONE truth test for every NT_* boolean env flag, in DSL.m because it loads first and Codegen.m
   shares this private context.

   Three mutually incompatible conventions used to coexist: presence-only (`=!= $Failed`, value
   ignored), exact-string (`=!= "1"`), and non-empty-string. So `NT_GEN_VERBOSE=0` turned verbosity
   ON, and — the one that actually mattered — `NT_NO_LABEL_CHECK=0` DISABLED the per-diagram label
   census, the guard that catches a Lorentz label occurring more than twice. That is not a knob: an
   uncaught repeated label is a silently wrong contraction, and the variable's own comment said
   "=1 disables", so setting it to 0 to be explicit did the opposite of what it read like.

   `SetEnvironment["VAR" -> None]` (the reset idiom the fixture generators use, e.g.
   gen_zaaqbq1_small_numeric.wls) makes Environment[] return $Failed, which must read as OFF — a
   naive `v =!= "0"` would call that truthy and leave the flag stuck ON for the rest of the session. *)
ntEnvFlag[name_String] :=
  With[{v = Environment[name]},
    StringQ[v] && MemberQ[{"1", "true", "yes", "on"}, ToLowerCase[StringTrim[v]]]];

(* ---- head classification ---------------------------------------------------- *)

(* A factor that participates in the tensor contraction (vs. a scalar coefficient).
   ntEpsilon (the Levi-Civita ε_{μνρσ}) is a Lorentz tensor produced by the γ5 trace theorem
   (gammaTraceSum5 in Codegen.m); it carries four Lorentz labels like a metric pair. *)
(* ntVec[q, mu] with a SYMBOLIC label is a tensor leg; ntVec[q, i_Integer] is the scalar component
   q_i (0-based, 0=temporal), resolved by the frame in the coeff like ntSP — so it is NOT a tensor.
   Express that by matching only a non-integer second arg here (an integer arg falls through to False). *)
tensorQ[_ntMetric | ntVec[_, Except[_Integer]] | _ntTransProj | _ntLongProj |
        _ntElectricProj | _ntMagneticProj | _ntSUNf | _ntSUNDeltaAdj |
        _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntSUNT | _ntSUNDeltaFund | _ntEpsilon |
        _ntSUNDiagFund | _ntSUNDiagAdj | _ntEpsFund | _ntDressedNum | _ntDiracSlot] = True;
tensorQ[_] = False;

(* The 4 SU(N) group heads carry their rank N as the FIRST argument, so colour SU(Nc),
   flavour/isospin SU(Nf), and any further group coexist in one network at the right rank.
   Lorentz, Dirac, and the SU(N) groups contract in disjoint index spaces (matched by shared
   label), so the label sets are treated uniformly while the builders stay distinct. *)
(* ntSUNDiag{Fund,Adj}[N, i, j, spec, scale] — a group δ that carries a PER-COMPONENT dressing.
   `spec` is a rules list {c -> name, ..., Default -> defName} of 1-based component indices to
   distinctly-named scalar dressing symbols (evaluated at the kinematic `scale`); unnamed components
   collapse to Default when given, else drop. They classify exactly like the matching plain δ
   (fundamental / adjoint) for index bookkeeping; the per-component dressing is folded numerically
   by sun_value_dressed at codegen time. *)
adjointSUNQ[_ntSUNf | _ntSUNDeltaAdj | _ntSUNDiagAdj] = True;  (* group-adjoint heads (bridge with Lorentz) *)
adjointSUNQ[_] = False;
fundamentalSUNQ[_ntSUNT | _ntSUNDeltaFund | _ntSUNDiagFund | _ntEpsFund] = True;   (* group-fundamental heads (quark-line) *)
fundamentalSUNQ[_] = False;

(* The SU(N) rank a group head builds against: its leading argument. *)
sunRankOf[h_] := First[h];

(* The contraction index labels carried by a tensor factor (momentum arg dropped). *)
labelsOf[ntMetric[mu_, nu_]]         := {mu, nu};
labelsOf[ntVec[_, mu_]]              := {mu};
labelsOf[ntTransProj[_, mu_, nu_]]   := {mu, nu};
labelsOf[ntLongProj[_, mu_, nu_]]    := {mu, nu};
labelsOf[ntElectricProj[_, mu_, nu_]] := {mu, nu};
labelsOf[ntMagneticProj[_, mu_, nu_]] := {mu, nu};
labelsOf[ntSUNf[_, a_, b_, c_]]      := {a, b, c};
labelsOf[ntSUNDeltaAdj[_, a_, b_]]   := {a, b};
labelsOf[ntGamma[mu_, din_, dout_]]  := {mu, din, dout};
labelsOf[ntGamma5[din_, dout_]]      := {din, dout};
(* ntSigma's OPEN Lorentz labels are its FREE legs (gluon ids); slashed legs carry a momentum, not an
   open id. The spinor axes din,dout are always open. *)
labelsOf[ntSigma[legA_, legB_, din_, dout_]] :=
  Join[Cases[{legA, legB}, {"free", mu_} :> mu], {din, dout}];
labelsOf[ntDeltaDirac[din_, dout_]]  := {din, dout};
labelsOf[ntSUNT[_, a_, i_, j_]]      := {a, i, j};
labelsOf[ntSUNDeltaFund[_, i_, j_]]  := {i, j};
labelsOf[ntSUNDiagFund[_, i_, j_, _]] := {i, j};  (* per-component-dressed fundamental δ *)
labelsOf[ntSUNDiagAdj[_, a_, b_, _]]  := {a, b};  (* per-component-dressed adjoint δ *)
labelsOf[ntEpsilon[a_, b_, c_, d_]]  := {a, b, c, d};
(* SU(N) fundamental Levi-Civita: N indices, all of them contraction labels. *)
labelsOf[h_ntEpsFund] := Rest[List @@ h];
(* ntDressedNum[options, din, dout] — a dressed propagator numerator kept EAGER (symbolic dressing
   collection). `options` is a tensor-FREE list of {coeffExpr, spec} (spec = {"ident"} | {"slash",vlc})
   so it exposes only its spinor in/out labels; the dressing collection is folded in C++ (one DPoly
   trace) instead of distributing the numerator into 2^D diagrams. *)
labelsOf[ntDressedNum[_, din_, dout_]] := {din, dout};
(* ntDiracSlot[options, din, dout, legs] — a collected Dirac slot (Stage 4 general form): a
   coefficient-weighted sum of Dirac structures sharing spinor in/out (din,dout) AND the SAME SET of
   open Lorentz legs `legs` (any count k>=0; k=0 = a propagator numerator, k>=1 = a vertex). `options`
   is a list of {coeffExpr, structureProduct}; the structure's free Lorentz ids are exactly `legs`, so
   the surrounding net closes them for every structure choice. It therefore exposes din,dout AND every
   open leg as real contraction labels. *)
labelsOf[ntDiracSlot[_, din_, dout_, legs_]] := Join[legs, {din, dout}];

(* Spinor (Dirac) axis labels of a head — the subset of labelsOf that lives in the
   spinor index space. Used to give spinor axes a disjoint id range so the et engine
   never contracts a spinor axis against a Lorentz/colour axis sharing an id. *)
spinorLabelsHead[ntGamma[_, din_, dout_]] := {din, dout};
spinorLabelsHead[ntGamma5[din_, dout_]]   := {din, dout};
spinorLabelsHead[ntSigma[_, _, din_, dout_]] := {din, dout};
spinorLabelsHead[ntDeltaDirac[din_, dout_]] := {din, dout};
spinorLabelsHead[ntDressedNum[_, din_, dout_]] := {din, dout};
spinorLabelsHead[ntDiracSlot[_, din_, dout_, _]] := {din, dout};
spinorLabelsHead[_]                       := {};
(* All spinor labels anywhere under a (sub)expression. *)
allSpinorLabels[e_] := DeleteDuplicates @ Flatten @ Cases[e, h_?tensorQ :> spinorLabelsHead[h], {0, Infinity}];

(* The momentum a factor carries (the source of runtime Var leaves), or None. *)
momentumOf[ntVec[q_, _]]             := q;
momentumOf[ntTransProj[q_, _, _]]    := q;
momentumOf[ntLongProj[q_, _, _]]     := q;
momentumOf[ntElectricProj[q_, _, _]] := q;
momentumOf[ntMagneticProj[q_, _, _]] := q;
momentumOf[_]                        := None;

(* Whether a momentum needs a 1/q^2 env slot (it sits inside a projector with a full denominator). *)
needsInvQ[ntTransProj[__]]    = True;
needsInvQ[ntLongProj[__]]     = True;
needsInvQ[ntElectricProj[__]] = True;  (* electric P_E = P_T − P_M uses 1/q² (and 1/|q⃗|², see needsInvSQ) *)
needsInvQ[_]                  = False;

(* Whether a momentum needs a SPATIAL 1/|q⃗|² env slot (finite-T electric/magnetic projectors). *)
needsInvSQ[ntElectricProj[__]] = True;
needsInvSQ[ntMagneticProj[__]] = True;
needsInvSQ[_]                  = False;

(* ---- self-trace normalization ----------------------------------------------- *)

(* A tensor with a repeated index is a self-trace (e.g. P^mu_mu). The et engine
   contracts pairwise BETWEEN tensors and never self-contracts one, so we relabel
   the second occurrence and insert the matching identity (Lorentz metric for a
   Lorentz index, adjoint delta for a colour index): P^mu_mu = P^{mu nu} d_{mu nu}.
   This keeps every index appearing on two distinct tensors, as the engine needs. *)
(* NOTE on Plus: tensorQ[Plus[...]] is False (head-matching with a catch-all), so a SUM vertex takes
   the `! tensorQ` branch and passes through untouched — a repeated index INSIDE a summand is then
   never split. That is safe only because a self-trace within one summand of an eager sum has never
   been produced by a flow; labelCensus's Plus branch would flag the resulting free-index mismatch if
   it were. Left as-is deliberately: relabelling inside a summand would have to keep every summand's
   free-index set aligned (the et::add precondition), which this factor-local rewrite cannot see. *)
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
scalarQ[e_] := FreeQ[e, _ntMetric | ntVec[_, Except[_Integer]] | _ntTransProj | _ntLongProj |
                       _ntElectricProj | _ntMagneticProj | _ntSUNf | _ntSUNDeltaAdj |
                       _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntSUNT | _ntSUNDeltaFund | _ntEpsilon |
                       _ntSUNDiagFund | _ntSUNDiagAdj | _ntEpsFund | _ntDressedNum | _ntDiracSlot];

(* The free (uncontracted) index labels of a tensor (sub)expression. A product sums
   indices that appear twice (free = appear once); a sum's summands share free indices
   (a vertex's legs); a bare head exposes all its labels. Used to group components and
   to align et::add operands. *)
freeIdx[e_] := Which[
  tensorQ[e],        labelsOf[e],
  Head[e] === Plus,  freeIdx[First[List @@ e]],
  (* A tensor^n is n copies sharing the SAME labels — a closed self-contraction (see compileTInv),
     so it exposes NO free index. Spelled out rather than left to the True branch below, which
     would return {} for the wrong reason and hide a malformed Power. *)
  Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 2 && ! scalarQ[e], {},
  (* count == 1 is free; count == 2 is contracted. checkLabels guarantees no label occurs more
     often, so the old OddQ test was equivalent — but OddQ would silently call a 4x label
     "contracted" and a 3x label "free" if the guard were ever bypassed. *)
  Head[e] === Times, Cases[Tally[Flatten[freeIdx /@ (List @@ e)]], {l_, c_} /; c == 1 :> l],
  True,              {}];

(* All index labels anywhere under a (sub)expression (free + internally summed) — the
   set that needs an axis id. *)
allLabels[e_] := DeleteDuplicates @ Flatten @ Cases[e, h_?tensorQ :> labelsOf[h], {0, Infinity}];

(* ---- connected components over shared free indices -------------------------- *)

(* Group tensor factors (heads OR Plus-vertices) into independent sub-networks: two
   factors are connected if their free indices intersect. Each component contracts to
   a scalar (colour stays a separate constant); keeps every contraction small. *)
connectedComponents[factors_List] := Module[{verts, edges, graph},
  verts = Range[Length[factors]];
  edges = Select[Subsets[verts, {2}],
    IntersectingQ[freeIdx[factors[[#[[1]]]]], freeIdx[factors[[#[[2]]]]]] &];
  graph = Graph[verts, UndirectedEdge @@@ edges];
  factors[[#]] & /@ (ConnectedComponents[graph])
];

(* Greedy contraction order: keep each successive factor sharing a free index with the
   running set, so contract_all's intermediates stay low-rank (never outer-product the
   whole thing). *)
orderFactors[fs_List] := Module[{rem = fs, out, used},
  If[fs === {}, Return[{}]];
  out = {First[rem]}; used = freeIdx[First[rem]]; rem = Rest[rem];
  While[rem =!= {},
    With[{pick = SelectFirst[rem, IntersectingQ[freeIdx[#], used] &, First[rem]]},
      AppendTo[out, pick]; used = Union[used, freeIdx[pick]];
      rem = DeleteCases[rem, pick, {1}, 1]]];
  out];

(* ---- sector-bridge expansion (keep colour and Lorentz contractions separate) ---- *)

(* Eager summation (et::add) is a win only when a structure-sum lives in ONE sector — a sum of
   Lorentz structures with a single colour factor pulled out (e.g. a 3-gluon vertex) stays one
   small Lorentz tensor. A SECTOR-BRIDGING sum is different: its colour sum is *correlated* with
   its Lorentz sum (each term pairs one colour structure with one Lorentz structure, e.g. a
   4-gluon vertex's f.f ⊗ metric.metric). et::add-ing those fuses the colour (dim N^2-1) and
   Lorentz (dim 4) axes into one ETensor whose entry count (and type) explodes — building it alone
   times out. The cure: DISTRIBUTE such a sum, so each term is colour/Lorentz index-disjoint again
   — colour folds to a constant, Lorentz stays small. That is linear in the (few) structures, not
   the product-of-sums blow-up. *)
sectorBridgeQ[p_] := (! FreeQ[p, _ntSUNf | _ntSUNDeltaAdj | _ntSUNDiagAdj]) &&
                     (! FreeQ[p, _ntMetric | ntVec[_, Except[_Integer]] | _ntTransProj | _ntLongProj |
                                 _ntElectricProj | _ntMagneticProj]);

(* The scalar (non-tensor) coefficient of a single summand. *)
scalarCoeffOf[t_] := Times @@ Select[If[Head[t] === Times, List @@ t, {t}], scalarQ];
(* A TENSOR-structure sum whose per-summand scalar coefficients are not all numeric must be
   DISTRIBUTED too: the et engine's eager add (et::add) can only fold NUMERIC per-structure
   scalars (it scales each structure by a compile-time Lit). A fermion propagator numerator
   Mq*deltaDirac + (.. dressings ..)*(gamma.vec) carries runtime dressing coefficients, so
   each Dirac structure must become its own diagram with that dressing as a scalar coeff.
   (Pure-scalar dressing sums — no tensor head — are left intact as coefficients.) *)
dressedStructureSumQ[p_Plus] := (! scalarQ[p]) && AnyTrue[List @@ p, ! NumericQ[scalarCoeffOf[#]] &];
(* Symbolic dressing collection: when $ntDressCollect is True a
   DIRAC dressed structure sum (a propagator numerator like Mq·δ + Z(p)·γ·p, all terms sharing the
   spinor in/out indices) is NOT distributed — it is kept eager and rewritten to a single ntDressedNum
   slot (rewriteDressedNums), so the diagram folds to ONE DPoly trace carrying its dressings as
   `dress` env leaves instead of exploding into 2^D diagrams. A colour↔Lorentz sectorBridge sum still
   distributes (sectors must split), and a non-Dirac dressed sum (no shared spinor indices) still
   distributes (the collection only handles the Dirac numerator case here). Default: $ntDressCollect
   False = today's full distribution. *)
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
(* A dressed Dirac sum is collectible only if every term is a propagator numerator we can ACTUALLY
   decompose (ident / slash, no open gluon leg). diracNumeratorSumQ is the cheap pre-filter (same 2
   open spinor indices, Dirac structure); dressedNumDecompose =!= $Failed is the exact test — a sum
   with an open Lorentz leg (a multi-structure quark-gluon / σ vertex) fails it and stays on the
   distribution path (so it never survives un-distributed AND un-rewritten → no eagernn abort). *)
(* Collectible if it is EITHER a propagator numerator (k=0, no open leg → ntDressedNum) OR a general
   Dirac vertex sum (k>=1 open legs → ntDiracSlot). Both require a NON-adjoint (quark-line) colour
   sector so colour factors out (sectorBridgeQ False) and the eager sum stays in one Lorentz/Dirac
   sector. The `=!= $Failed` decompose is the exact test in each case. *)
(* VERTEX collection (k>=1 open-leg ntDiracSlot) is OPT-IN — OFF by default. On a high-combination
   flow (full-basis ZAAqbq: ~10^6 structure x dressing combinations per net) the codegen expansion
   materialises the Cartesian product and OOMs the WolframKernel, so such flows MUST distribute (the
   pre-plan behaviour). The pre-existing PROPAGATOR collection (k=0) is bounded (2 options/propagator)
   and stays ON unconditionally. Enable the vertex path with NT_VERTEX_COLLECT=1 (or set
   $ntVertexCollect=True) on the small-P flows where it wins (e.g. za3_147: 8.9 s vs 13.3 s). *)
$ntVertexCollect = ntEnvFlag["NT_VERTEX_COLLECT"];
(* ALL-NUMERIC Dirac sums as slots (opt-in). `dressedStructureSumQ` demands at least one NON-numeric
   summand coefficient, because it was written to decide whether a DRESSING sum must be distributed.
   A multi-term PROJECTOR is collateral damage: its coefficients are pure numbers, so it fails that
   test, is therefore never collectible, and survives as a raw Plus all the way to
   splitColourGroupsInv — where `Expand[Times @@ needExpand]` materialises its Cartesian product with
   the rest of the diagram, once per diagram, per generation.
   That is measurably the dominant cost of a projector-heavy flow: transAAqbqMinimal's element 2
   (~10 Dirac terms, 30 gammas) against element 1 (3 terms, 4 gammas) is 1056 vs 14 distinct emitted
   DiracNets and 328 s vs 62 s of generation, from two otherwise identical flow definitions.
   Such a sum satisfies every OTHER condition of the slot path — shared open spinor pair, shared open
   Lorentz set, colour factoring out (it sits outside the Plus) — so relaxing the gate for the SLOT
   disjunct alone turns it into one chain with an N-option slot.
   Scoped deliberately: the k=0 propagator-numerator disjunct keeps `dressedStructureSumQ`, so an
   all-numeric propagator numerator still distributes exactly as today. And `distributeQ` below is
   untouched by construction — it tests `dressedStructureSumQ` first, which is False for these sums,
   so its verdict (do not distribute) is the same either way. *)
(* DEFAULT ON, and deliberately NOT gated behind $ntVertexCollect. The two are different risks that
   happen to share a mechanism. $ntVertexCollect is opt-in because a DRESSED vertex sum multiplies the
   structure count by the dressing count — the ~10^6 structure x dressing combinations per net that
   OOMed the kernel. An all-numeric projector sum has no dressing dimension at all (its options carry
   no dress atoms), so its combination count is bounded by the projector's own term count: measured on
   ZAAqbq2 the sub-term count is IDENTICAL either way (12,721,032) — the work moves, it does not grow.
   Hatch names the legacy path, per the convention NT_DIRAC_FLAT set. *)
$ntSlotCollectNumeric = ! ntEnvFlag["NT_NO_SLOT_COLLECT_NUMERIC"];
collectibleDiracSumQ[p_Plus] := ! sectorBridgeQ[p] &&
  ((dressedStructureSumQ[p] && diracNumeratorSumQ[p] && dressedNumDecompose[p] =!= $Failed) ||
   ((TrueQ[$ntVertexCollect] && dressedStructureSumQ[p]) ||
      (TrueQ[$ntSlotCollectNumeric] && ! dressedStructureSumQ[p])) &&
     diracSlotSumQ[p] && diracSlotDecompose[p] =!= $Failed);
collectibleDiracSumQ[_] := False;
distributeQ[p_] := sectorBridgeQ[p] ||
  (dressedStructureSumQ[p] && ! (TrueQ[$ntDressCollect] && collectibleDiracSumQ[p]));

(* The SET of γ-count parities the diagram's terms would carry IF every EAGER Dirac Plus were
   distributed — {0}, {1}, or {0,1}. Used for the odd-trace-vanishing verdict: a closed trace is
   identically zero only when EVERY branch has an odd γ count, i.e. the set is exactly {1}.
   A diagram-global `Count[_ntGamma, Infinity]` cannot decide this — it sums γ's ACROSS the branches
   of an eager Plus and so reports a parity no actual term has (e.g. a projector 1 + p̸ against two
   γ's counts 3 → "odd" → the whole diagram is dropped, taking its non-vanishing even branch with it).
   Computed WITHOUT Expand: the sets are capped at two elements and fold pairwise, so this is linear
   in the factor count rather than a product-of-sums blow-up.
   Only _ntGamma is counted, matching the old test: ntSigma (2 γ) and ntDeltaDirac (0 γ) are both
   EVEN and so correctly contribute parity 0. ntDressedNum cannot appear yet (rewriteDressedNums runs
   later, inside analyseDiagram). *)
diracParities[e_Plus] := Union @@ (diracParities /@ (List @@ e));
diracParities[e_Times] := Fold[Union[Flatten[Mod[Outer[Plus, #1, #2], 2]]] &,
                               {0}, diracParities /@ (List @@ e)];
(* b^n is n copies of b multiplied. Do NOT write this as diracParities[Times @@ ConstantArray[b, n]]:
   Times immediately re-collapses the copies back to b^n, so that recurses until $RecursionLimit and
   returns garbage — and it fires on ORDINARY SCALAR powers (sp[p2,p3]^2 in a Gram determinant), not
   just tensor powers. Fold the parity arithmetically instead: n copies of a definite parity p give
   n*p mod 2, and if b itself has both parities then so does any power of it. *)
diracParities[Power[b_, n_Integer?Positive]] := With[{s = diracParities[b]},
  If[Length[s] > 1, {0, 1}, {Mod[n * First[s], 2]}]];
(* {0, Infinity}, NOT Infinity: this catch-all is reached with a BARE factor (a single ntGamma head
   pulled out of the enclosing Times), and level spec Infinity means {1, Infinity} — it would skip
   level 0 and count that γ as zero. The old diagram-global test never hit this because it always
   ran on a whole Times, where every γ sits at level 1. *)
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

dressedNumDecompose[p_Plus] := Module[
   (* Flatten only this local numerator sum. At finite T, psdash[p] contains
      gamma.mu vecs[p, mu]; after fixed-component normalization the temporal
      subtraction must become a separate slash option rather than making the
      complete propagator numerator non-collectible. *)
   {rows = dressedNumTerm /@ (List @@ Expand[p]), din, dout, others, common, scalFacs, commonScal, opts},
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

(* Inverse of the ntDressedNum rewrite: expand a collected numerator back into its distributed Dirac
   structure sum (Σ coeff_i · {ntDeltaDirac | ntGamma·Σc·ntVec}). Used when re-distributing a
   collected diagram (the small-D gate in NumTrace). Each slash gets a fresh Lorentz label for its
   contracted leg. *)
expandDressedNum[ntDressedNum[opts_, din_, dout_]] := Plus @@ (Function[opt,
  opt[[1]] * Switch[opt[[2, 1]],
    "ident", ntDeltaDirac[din, dout],
    "slash", With[{mu = Unique["dexp"]},
       ntGamma[mu, din, dout] * (Plus @@ ((#[[1]] * ntVec[#[[2]], mu]) & /@ opt[[2, 2]]))]]] /@ opts);

(* ---- general collected Dirac slot (Stage 4, ANY open-leg count) --------------------------------
   A collected Dirac slot is a coefficient-weighted sum of Dirac structures that all share the spinor
   in/out pair AND the SAME SET of open Lorentz legs `{μ...}` (k>=0). k=0 is the propagator numerator
   (handled by the ntDressedNum path); k>=1 is a vertex with k open gluon legs (Aqbq: 1; AAqbq: 2; …).
   The C++ engine already closes an arbitrary number of open legs against the net, so the front end
   only has to keep the sum EAGER and package each structure as one option. Unlike the propagator
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
diracSlotSumQ[p_Plus] := Module[{terms = List @@ Expand[p], opens, lors},
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
diracSlotDecompose[p_Plus] := Module[
  {terms = List @@ Expand[p], legs, opens, din, dout, io, ins, outs, rows, cols, common, scals, commonScal, opts},
  If[! diracSlotSumQ[p], Return[$Failed]];
  legs   = Sort @ openLorentzOf[First[terms]];
  opens  = openSpinorOf[First[terms]];
  (* ORIENTED din/dout: the chain runs din→dout. Each Dirac head is an (in,out) spinor edge; the open
     in-leg is the open spinor label that is some head's `in` but no head's `out` (dout is the reverse).
     Orientation matters — the option's tokens are spliced in chain order, so a reversed din/dout would
     emit the trace backwards. *)
  io = Cases[If[Head[First[terms]] === Times, List @@ First[terms], {First[terms]}],
        ntGamma[_, a_, b_] | ntGamma5[a_, b_] | ntDeltaDirac[a_, b_] | ntSigma[_, _, a_, b_] :> {a, b}];
  ins = io[[All, 1]]; outs = io[[All, 2]];
  din  = SelectFirst[opens, MemberQ[ins, #] && ! MemberQ[outs, #] &, First[opens]];
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

(* Inverse: expand a collected slot back into its distributed Dirac structure sum (for redistDiagram's
   small-D cross-check). Each option is coeff × its whole structure product, so this is exact. *)
expandDiracSlot[ntDiracSlot[opts_, _, _, _]] := Plus @@ ((#[[1]] * #[[2]]) & /@ opts);

(* Re-distribute ONE analysed (collected) diagram back to the non-collected path: rebuild its net,
   expand every ntDressedNum into the Dirac structure sum, distribute, drop odd-gamma traces, and
   re-analyse (with collection OFF so no ntDressedNum is recreated). Used by the small-D gate in
   NumTrace. *)
redistDiagram[diag_] := Block[{$ntDressCollect = False},
  Module[{net = diag["Coeff"] * Times @@ Flatten[(#["Factors"] &) /@ diag["Components"]]},
    net = net /. nd_ntDressedNum :> expandDressedNum[nd];
    net = net /. ds_ntDiracSlot :> expandDiracSlot[ds];
    With[{ex = expandBridges[net]},
      analyseDiagram /@ Select[If[Head[ex] === Plus, List @@ ex, {ex}],
        ! vanishingOddTraceQ[#] &]]]];

(* Distribute every colour<->Lorentz-bridging sum into its surrounding product (only that
   sum; single-sector sums are left intact for et::add). Turns a bridging diagram into a
   small linear sum of sector-separable diagrams. Explicit recursion (not a //. rule over
   Orderless Times, which backtracks catastrophically on a large net). *)
expandBridges[e_Plus] := Plus @@ (expandBridges /@ (List @@ e));
expandBridges[e_Times] := Module[{factors = List @@ e, bridge},
  bridge = FirstCase[factors, p_Plus /; distributeQ[p], Missing[]];
  If[MissingQ[bridge],
    Times @@ (expandBridges /@ factors),
    expandBridges[Plus @@ (Times @@ Append[DeleteCases[factors, bridge, {1}, 1], #] & /@ (List @@ bridge))]]];
(* Times is Flat+Orderless and auto-collects identical factors, so two byte-identical bridging sums
   become Power[sum, 2] BEFORE we ever see them. The selector above matches p_Plus only, so such a
   Power is never chosen as a bridge and would fall through the catch-all below UN-DISTRIBUTED —
   fusing the colour and Lorentz axes into one giant ETensor (the exact blow-up expandBridges
   exists to prevent), or worse being read by compileTInv as a closed self-contraction. It cannot
   simply be re-expanded here: Times would immediately re-collapse the copies unless their labels
   were freshened first, and two vertices that legitimately share every label are themselves a bug.
   No known flow produces this; fail loudly if one ever does. *)
NumTrace::bridgepow = "expandBridges: a colour<->Lorentz-bridging sum appears raised to the power \
`1`, i.e. as `1` byte-identical factors sharing every index label. Such a sum cannot be distributed \
(and identical labels on two distinct vertices are themselves malformed). Offending base:\n`2`";
expandBridges[e : Power[b_Plus, n_Integer]] /; n >= 2 && distributeQ[b] :=
  (Message[NumTrace::bridgepow, n, Short[b, 6]]; Abort[]);
expandBridges[e_] := e;

(* ---- fixed Lorentz components (the finite-T γ0/γi split) ---------------------

   The four-quark Fierz bases (and any finite-T 3+1 split) are written with Lorentz indices pinned
   to a CONCRETE component: ntGamma[0, d1, d2] is γ^0, not "γ^μ with a label named 0". The DSL's
   label machinery has no such notion — labelsOf[ntGamma[mu,...]] returns {mu,...} unconditionally
   — so that literal 0 was read as an ordinary contraction label, appeared in a dozen factors at
   once, and checkLabels aborted with the (correct but baffling) `privclash`.

   The fix is a REWRITE, not new machinery: a fixed component is a contraction with the constant
   unit basis vector e_i, so

       γ^i          ->  ntGamma[μ, d1, d2] ntVec[ntUnitVec[i], μ]     (μ a fresh private dummy)
       g^{i ν}      ->  ntVec[ntUnitVec[i], ν]
       g^{i j}      ->  δ_ij                                          (Euclidean metric)
       σ's free leg ->  the existing SLASH leg against ntUnitVec[i]

   ntUnitVec[i] is an ordinary MOMENTUM symbol (momentumOf picks it up, buildEnv gives it a Base),
   whose frame components NumTrace injects as UnitVector[4, i+1]. Downstream nothing changes:
   compileDirac's vecOf already turns a γ whose μ carries an ntVec into `dslash({{1.0, base}})`,
   which is exactly γ contracted with that vector. So this needs NO new DFac kind and NO C++ change
   — a fixed-component γ IS a slash, and frameMask prunes the three zero components for free.

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

   vecs[q, mu] is the spatial part of q as a 4-vector: components {0, q_1, q_2, q_3}. FromFunKit
   rewrites it to ntVec[ntSpatialVec[q], mu], so a SPATIAL SLASH vecs[q,mu] gamma[mu,d1,d2] is an
   ordinary slash against a different momentum — compileDirac's vecOf emits `dslash({{1.0, base}})`
   for it unchanged. Exactly the ntUnitVec trick above: a new momentum LEAF, no new DFac kind, no C++
   change. The zero temporal component costs nothing either: it is a structural zero in the frame
   spec's component table, so it never becomes a term. (It also lowers the momentum's frameMask, but
   do not read anything into that — the numeric backend's `lvec<Lbl,Base,Mask>` ignores its Mask
   argument entirely, and a slash carries no mask at all.)

   Contrast ntSPS, the spatial scalar PRODUCT: that is a scalar coefficient the frame resolves as a
   components-1..3 dot, so it needs no leaf at all. The two must agree, which the codegen gate
   gen_spatialvec_numeric.wls pins.

   LINEARITY first. The spatial projection is linear, so push it through sums and numeric factors
   (FormTracer does the same, FormTracer.m:888-889) BEFORE anything else looks at the momentum. Two
   reasons: (i) only BASE momenta then become frame keys, so ntSpatialVec[p - l] does not mint a
   third leaf whose components duplicate those of ntSpatialVec[p] and ntSpatialVec[l]; (ii) the leaf
   is then an atom by the time canonicalizeMomentumSigns' negMomQ inspects it, so the sign
   convention applies to the spatial momentum itself rather than to whatever sat inside it. *)
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

   These entries make the spatial vector a first-class frame/env citizen — buildEnv gives it a Base,
   resolveComponents resolves it, frameMask masks it. They are NOT what the numeric backend computes
   with (measured by mutation: perturbing them leaves every emitted kernel byte-identical, because
   `lvec` drops its Mask and a slash never had one): the component table comes from the frame SPEC
   (Codegen.m unitLoopFrameSpec /
   unitLoopMixedFrameSpec / polyFrameSpec), where a spatial vector is likewise derived from its
   parent with slot 1 zeroed — sharing the parent's ntU$ unit group rather than minting a duplicate.
   Two derivations, one invariant ("the parent's components, temporal slot zeroed"), which is why
   they cannot drift apart. *)
spatialVecFrame[net_, frame_] := Association[
  Function[sv, Module[{c = resolveComponents[First[sv], frame]},
      If[! MatchQ[c, {_, _, _, _}],
        Message[NumTrace::spatialframe, First[sv], c]; Abort[]];
      sv -> ReplacePart[c, 1 -> 0]]] /@
    DeleteDuplicates[Cases[net, _ntSpatialVec, {0, Infinity}]]];

(* ---- SU(N) FUNDAMENTAL Levi-Civita -------------------------------------------

   ntEpsFund[N, i1, ..., iN] is the totally antisymmetric invariant of SU(N) in the FUNDAMENTAL
   space, so it carries exactly N indices (colour SU(3): 3; isospin SU(2): 2). It reaches us from
   FunKit's epsFundCol/epsFundFlav, which the four-quark Fierz bases use for their diquark /
   colour-antisymmetric channels.

   Like the fixed-component gamma, this is a REWRITE into primitives the engine already has, not a
   new engine token: an epsilon is contracted ONLY in pairs (a lone epsilon is not an SU(N)
   invariant), and a pair folds to a determinant of Kronecker deltas,

       eps_{a1..ak c1..cm} eps_{a1..ak d1..dm}  =  k! * det( delta_{c_p d_q} )        (m = N - k)

   which is exactly ntSUNDeltaFund. So there is NO C++ change: SUNFac's fixed 3-slot layout never
   has to represent an N-index object. The result is a Plus of delta products with numeric
   coefficients, which lands on the constant-colour branch-list path (compileColGSum in Codegen.m).

   Applied BEFORE expandBridges/checkLabels so the object those validate is the one that compiles.

   WHY THE HEAD IS DELIBERATELY NOT REGISTERED IN Codegen's colour tables (ctHeadsInv / colFacG /
   labelDimAssoc): an ntEpsFund that somehow survives this rewrite must FAIL, not be emitted. It then
   trips colFacG's catch-all (MakeNTKernel::colleak) and, behind that, ntExportCpp's nt*-head regex. *)

(* The dimension of the index space an epsilon head lives in. One line today; the single point a
   hypothetical adjoint epsilon (dimension N^2-1) would extend. *)
epsDimOf[ntEpsFund[n_, __]] := n;

(* The pair contraction, written DIMENSION-PARAMETRIC (dim and the delta constructor are the only
   things that were N-specific). That makes it unit-testable at dim = 2..5 against LeviCivitaTensor
   without any head existing at those dimensions — see the gate's brute-force oracle. *)
$ntEpsMaxPairTerms = 720;   (* 6!; 7! = 5040 would breach $ntColSumMaxBranches (4096) on its own *)
NumTrace::epsbig = "expandFundEps: an epsilon pair in dimension `1` sharing `2` index/indices \
expands to `3`! = `4` Kronecker-delta terms (limit `5`). Emitting these would hand the colour \
branch-list lowering a list it would either reject far downstream with an opaque branch count, or \
— worse — accept and turn into a generator source large enough to OOM the C++ compiler, with \
nothing pointing at a Levi-Civita as the cause. Contract more indices between the two epsilons, or \
raise $ntEpsMaxPairTerms if a flow genuinely needs this.";

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
expandFundEps[e_] := Module[{res},
  res = expandFundEpsRec[e];
  With[{left = Cases[res, _ntEpsFund, {0, Infinity}]},
    If[left =!= {}, Message[NumTrace::epsodd, Length[left], Short[left, 4]]; Abort[]]];
  res];

expandFundEpsRec[e_Plus] := expandFundEpsRec /@ e;
(* eps^2 is a self-contraction of an epsilon with itself: every index is shared, so k = N, m = 0 and
   the value is N!. Compute it ARITHMETICALLY — do NOT write this as expandFundEpsRec[b b], because
   Times immediately re-collapses two identical factors back to b^2 and the rule recurses until
   $IterationLimit. (Exactly the trap the diracParities Power clause documents.) Anything above the
   square is an odd/malformed count and falls through to the epsodd guard. *)
expandFundEpsRec[Power[b_ntEpsFund, 2]] := Module[{idx = Rest[List @@ b], d = epsDimOf[b]},
  If[Length[idx] =!= d, Message[NumTrace::epsrank, Length[idx], d, b]; Abort[]];
  If[Length[DeleteDuplicates[idx]] =!= Length[idx], 0,
     epsPairExpand[d, idx, idx, ntSUNDeltaFund[d, #1, #2] &]]];
expandFundEpsRec[e_Times] := Module[{fs, eps, rest, cand, pair, uu, vv, plusEps, host, others},
  (* recurse into the factors FIRST: an epsilon pair frequently lives inside a Plus factor (a
     multi-term projector), and pairing only at this level would leave it untouched — the survivors
     then trip the epsodd guard, which is correct but unhelpful. *)
  fs = expandFundEpsRec /@ (List @@ e);
  eps = Cases[fs, _ntEpsFund];
  (* Times @@ fs, NOT e: the factors were just recursed into (a Plus factor typically had its own
     pairs contracted), so returning the original here would silently discard that work. The test is
     FreeQ over all of fs, not just the bare factors: an epsilon may still sit inside a Plus factor
     awaiting a partner from out here (the straddle case handled at the bottom). *)
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
     Candidates are restricted to EQUAL RANK. Pairing two epsilons of different rank is not a weaker
     identity, it is nonsense: epsPairExpand would build a determinant sized by one partner's rank
     and silently DROP the surplus indices of the other. That is exactly how a rank-2 flavour epsilon
     got contracted against a rank-3 colour one in the four-quark Fierz diquark vertex — one colour
     leg vanished, the two summands of the eager vertex sum ended up exposing DIFFERENT free indices,
     and the failure surfaced far downstream as NumTrace::plusfree. Unequal ranks are left unpaired
     here and reported by the epsodd guard, which names them. *)
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
       rank left — under a well-formed input each space carries an even count, so two survivors of
       one rank must be partners). Otherwise refuse. *)
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
     the diquark vertex of the four-quark Fierz bases. eps.eps = k! det(delta) is only applicable
     where both partners multiply each other, so distribute the remaining factors into ONE such sum
     and recurse; each step removes one epsilon-bearing Plus, so this terminates.
     Times, never Expand: a summand's own internal sums (the Dirac structure sum here) stay intact,
     so this costs the epsilon-bearing sum's width and nothing more. *)
  If[eps =!= {} || ! FreeQ[rest, _ntEpsFund],
    plusEps = Select[rest, Head[#] === Plus && ! FreeQ[#, _ntEpsFund] &];
    If[plusEps =!= {},
      host   = First[plusEps];
      others = Times @@ Join[DeleteCases[rest, host, {1}, 1], eps];
      Return[Plus @@ (expandFundEpsRec[others #] & /@ (List @@ host)), Module]]];
  Times @@ Join[rest, eps]];   (* a leftover odd epsilon rides along to the epsodd guard *)
expandFundEpsRec[e_] := e;

(* ---- momentum sign canonicalisation ----------------------------------------- *)

(* buildEnv keys momentum Bases on the momentum EXPRESSION, so `q` and `-q` used to get separate
   Bases and separate Inv slots — measured on with_mesons' lambda1L3D: 6 exact ± pairs among 16
   bases, and 14 inv slots holding only 7 distinct values. Worse than the wasted slots, the sign
   twins block the global sub-term dedup and the emitted-body dedup (both key on the emitted net
   source, which names the base): the two-momentum flow deduped only 7% of its trace bodies vs 30%
   for its fast-path sibling. Every momentum-carrying head is either LINEAR in its momentum (a
   vector/slash leg: vec(-q) = -vec(q), the -1 surfacing as an ordinary scalar factor in the
   diagram's Times/Plus tree) or EVEN (the projectors and their 1/q² inv atoms), so a canonical
   sign representative is exact. The convention: flip iff the coefficient of the FIRST variable
   (canonical Sort order) is numerically negative — deterministic, and maps each ± pair to one rep.

   NT_NO_SIGN_CANON: escape hatch back to expression-keyed bases; exists for grading (the
   lambda3d_small control kernel) and as the rollback. *)

negMomQ[q_] := Module[{vars = Sort[Variables[q]], c},
  vars =!= {} && (c = Coefficient[q, First[vars]]; NumericQ[c] && c < 0)];

canonSlashPairs[vlc_List] :=
  Replace[vlc, {c_, q_} /; negMomQ[q] :> {-c, Expand[-q]}, {1}];

canonicalizeMomentumSigns[net_] :=
  If[ntEnvFlag["NT_NO_SIGN_CANON"],
    net,
    net /. {
      ntVec[q_, l_] /; negMomQ[q] :> -ntVec[Expand[-q], l],
      ntTransProj[q_, a_, b_] /; negMomQ[q] :> ntTransProj[Expand[-q], a, b],
      ntLongProj[q_, a_, b_] /; negMomQ[q] :> ntLongProj[Expand[-q], a, b],
      ntElectricProj[q_, a_, b_] /; negMomQ[q] :> ntElectricProj[Expand[-q], a, b],
      ntMagneticProj[q_, a_, b_] /; negMomQ[q] :> ntMagneticProj[Expand[-q], a, b],
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
  {env, inv}  (* inv is now the total env size NEnv *)
];

(* Per-momentum component mask: bit i set <=> component i is structurally nonzero in
   the frame. Drives the et builder Mask template arg (zero components prune away). *)
frameMask[components_List] := FromDigits[Reverse[Boole[# =!= 0 && # =!= 0.] & /@ components], 2];

(* ---- NumTrace --------------------------------------------------------------- *)

(* NOTE (do not re-attempt): running analyseDiagram / labelCensus over Wolfram SUBKERNELS was tried
   and REMOVED. Both are pure per-diagram maps, so parallelism is trivially correct — but it is
   marshalling-bound, not compute-bound: the kernel serialises every diagram AND its result across
   processes, costing about as much as the work. Measured on full-basis ZAAqbq (3350 diagrams):
   NumTrace 55 s serial vs 55 s parallel, and the label census got *slower* (it is too cheap to beat
   the transfer). It also had a correctness trap — the subkernels needed every collection global
   ($ntDressCollect, $ntVertexCollect, …) pushed to them by hand, and a flow that sets one by symbol
   assignment rather than by environment variable would silently analyse differently on the workers.
   The real front-end lever is NT_NO_LABEL_CHECK (skip the census, ~14%). See
   NUMTRACER_TRACE_PERF_FINDINGS.md. *)

(* Each SU(N) group head carries its own rank N as the first argument (baked in when the
   network is built — Global`Nc for colour, the FromFunKit "FlavourGroup" option for the
   isospin group), so NumTrace itself takes no group option. *)
Options[NumTrace] = {"Frame" -> <||>, "Args" -> {}, "Dressings" -> {}, "DressingCollection" -> True};

NumTrace[net_, OptionsPattern[]] := Module[
  {frame, args, dress, badRanks, net2, diagrams, allMom, invMom, invSMom, env, nenv, diags, ntT0},
(* WHOLE-NumTrace wall clock. The stage timers below ([prof] expandBridges / checkLabels /
   analyseDiagram) do not add up to the call: canonicalizeMomentumSigns, expandFundEps,
   expandSpatialVecs, expandFixedComponents, the odd-trace prune and buildEnv are all untimed.
   Without the total there is no way to tell "the untimed remainder is noise" from "the untimed
   remainder IS the flow's cost", which is exactly the question a profiling run is asked. *)
  ntT0   = AbsoluteTime[];
  frame  = OptionValue["Frame"];
  args   = OptionValue["Args"];
  dress  = OptionValue["Dressings"];
  (* symbolic dressing collection: keep dressed Dirac numerators eager (one DPoly trace) instead of
     distributing into 2^D diagrams. Default False = today's behaviour. Set here so expandBridges
     (distributeQ) and analyseDiagram (rewriteDressedNums) both see it; each NumTrace call sets it
     from its option (default False), so it never leaks across kernels. *)
  $ntDressCollect = TrueQ[OptionValue["DressingCollection"]];

  (* SU(N) ranks must be compile-time integers (they pick the correct-dimension typed-out
     group matrices). Every group head's leading argument N is checked up front, so an
     undefined symbol aborts with a clear message rather than being baked into the kernel. *)
  badRanks = DeleteDuplicates @ Select[
    Cases[net, h : (_ntSUNf | _ntSUNDeltaAdj | _ntSUNT | _ntSUNDeltaFund | _ntSUNDiagFund | _ntSUNDiagAdj | _ntEpsFund) :> sunRankOf[h], Infinity],
    ! (IntegerQ[#] && # >= 1) &];
  If[badRanks =!= {},
    Print["NumTrace: every SU(N) head must carry an integer rank N >= 1 as its first argument ",
          "(got ", badRanks, "). Set Global`Nc (SetNc[3]) before tracing colour heads, and pass ",
          "\"FlavourGroup\" -> n to FromFunKit for the isospin group."]; Abort[]];
  (* A closed quark flavour loop folds to Global`Nf. When SetNf[n] has been called Global`Nf is
     an integer and has already folded into the coefficients (nothing to check). Only when it is
     still the unset symbol AND it appears in the network is the flavour count genuinely undefined
     — check for the symbol then (guarding so the integer value is not mistaken for a literal). *)
  If[! IntegerQ[Global`Nf] && ! FreeQ[net, Global`Nf],
    Print["NumTrace: flavour Nf is not a defined integer (symbolic Nf present in the network). ",
          "Call SetNf[2] (TensorBases) before generating."]; Abort[]];

  (* FIXED LORENTZ COMPONENTS (γ^0 & co, the finite-T 3+1 split used by the four-quark Fierz bases)
     are rewritten into contractions with constant unit basis vectors FIRST, so that no integer
     Lorentz slot ever reaches the label machinery below. See expandFixedComponents. The unit
     vectors join the frame as ordinary momenta, so every existing frame builder stays untouched. *)
  net2 = canonicalizeMomentumSigns @ expandFundEps @ expandSpatialVecs @ expandFixedComponents[net];
  With[{bad = DeleteDuplicates @ Cases[net2, ntUnitVec[i_] :> i, {0, Infinity}]},
    If[! AllTrue[bad, IntegerQ[#] && 0 <= # <= 3 &],
      Message[NumTrace::fixcomp, Select[bad, ! (IntegerQ[#] && 0 <= # <= 3) &]]; Abort[]]];
  frame = Join[frame, unitVecFrame[net2]];
  (* FINITE-T SPATIAL VECTORS (FormTracer's vecs) join the frame the same way, but AFTER the unit
     vectors — a spatial vector's components are read off its parent's, so the parent must already
     resolve. See expandSpatialVecs / spatialVecFrame. *)
  frame = Join[frame, spatialVecFrame[net2, frame]];

  (* the top-level sum is the (linear) sum of DIAGRAMS; keep single-sector vertex sums eager
     via et::add, but distribute colour<->Lorentz-bridging sums so the two sectors never fuse
     into one giant ETensor. *)
  ntLog["[prof] NumTrace expandBridges: ", First@AbsoluteTiming[
  diagrams = With[{ex = expandBridges[net2]}, If[Head[ex] === Plus, List @@ ex, {ex}]];], " s"];
  (* odd-trace vanishing: a closed Dirac trace of an ODD number of gammas is identically zero
     (and would otherwise carry a spurious imaginary I^odd coefficient from the vertex/basis
     normalizations). Drop such diagrams. (gamma5-bearing traces are exempt from the rule.)
     The verdict is PER BRANCH (diracParities), not on a diagram-global gamma count: an EAGER Dirac
     Plus left standing by expandBridges can carry a VARIABLE gamma parity, and summing the count
     across its branches yields a parity no actual term has. Two ways that bites:
       * a COLLECTIBLE dressed propagator sum (Mq·δ + Z·γ·p̸, kept eager under collection) — its δ
         ("ident") branch has 0 gammas, its γ·p̸ ("slash") branch 1. (e.g. struct-4 `l̸1 γ^ρ` + a
         propagator counts 3 γ → "odd" → the diagram was dropped, killing its non-zero Mq·δ branch
         — the 1/4/7 collection bug.)
       * a MULTI-TERM PROJECTOR's Plus — all-numeric coefficients and single-sector, so distributeQ
         is False (left eager) AND collectibleDiracSumQ is False (that predicate needs a NON-numeric
         coefficient). It was exempted by neither, so the oblique-metric dual projector for e.g.
         AqbqDirect8 structure 7 traced to an identically-zero kernel with no error raised.
     Only an ALL-odd-branch diagram is dropped; a mixed-parity one is KEPT, and the per-branch
     matrix-product trace zeroes its odd-gamma branches on its own (splitColourGroupsInv expands any
     γ-bearing Plus and compiles each branch separately, so this filter is purely a pruning
     optimisation sitting upstream of already-correct code). *)
  diagrams = Select[diagrams, ! vanishingOddTraceQ[#] &];

  (* Validate every distributed diagram BEFORE analyseDiagram assigns axis ids: it gives ONE id
     per DISTINCT label (see below), so a label occurring 4x becomes four axes sharing an id and
     the et engine mis-pairs them into a silently wrong number. Checked per diagram (not on the
     raw net): only after expandBridges is each diagram a flat Times in which "1 = free,
     2 = contracted" is the actual invariant. *)
  ntLog["[prof] NumTrace checkLabels: ", First@AbsoluteTiming[
  With[{frees = If[TrueQ[$ntCheckLabels],
      (* the census (labelCensus) is pure; the abort/Message validation is kept separate so a
         failure reports the diagram index (see checkLabels) *)
      With[{census = labelCensus /@ diagrams},
        MapThread[checkLabels[#1, #2, #3] &, {diagrams, census, Range[Length[diagrams]]}]],
      ConstantArray[{}, Length[diagrams]]]},
    ntLog["[labels] ", Length[diagrams], " diagram(s) validated; free-index set(s) = ",
      DeleteDuplicates[Sort /@ frees]]];], " s"];

  (* global env layout: every distinct momentum, and which ones need a 1/q^2 slot *)
  (* net2, not net: the unit basis vectors introduced by expandFixedComponents — and the spatial
     vectors introduced by expandSpatialVecs — are ordinary momenta and MUST get an env Base, or
     compileDirac's slash emission finds them absent. *)
  allMom = DeleteDuplicates @ Cases[net2, f_?tensorQ :> momentumOf[f], Infinity] // DeleteCases[None];
  invMom = DeleteDuplicates @ Cases[net2, f_?(needsInvQ) :> momentumOf[f], Infinity];
  invSMom = DeleteDuplicates @ Cases[net2, f_?(needsInvSQ) :> momentumOf[f], Infinity];
  {env, nenv} = buildEnv[allMom, invMom, invSMom];

  ntLog["[prof] NumTrace analyseDiagram (", Length[diagrams], " diagrams): ",
    First@AbsoluteTiming[diags = analyseDiagram /@ diagrams], " s"];

  (* ---- NO FLAVOUR DELTA MAY LEAVE HERE -------------------------------------------------------
     This is the first point at which the fundamental-flavour residue is final: both contractFlavour
     passes have run and promoteFlavResidue has had its chance. A flavDelta that is still standing is
     neither contracted nor in the engine, and it fails SILENTLY downstream — scalarQ is a FreeQ over
     the nt* heads, so it is classified as a scalar COEFFICIENT, its indices become invisible to
     labelsOf/freeIdx/checkLabels, and CForm prints it into the kernel as
     `NumTracer_Private_flavDelta(F1, F2)`: a bare identifier both GCC and Clang accept as an
     undeclared call. Refuse it here, where the diagram is still nameable, rather than three layers
     away in a compiler error. (The textual $ntCppLeakPatterns entry in Codegen.m is the backstop for
     the same class; this one gives the better message.) *)
  With[{leak = DeleteDuplicates @ Cases[diags, _flavDelta, {0, Infinity}]},
    If[leak =!= {}, Message[NumTrace::flavleak, Short[leak, 8]]; Abort[]]];

  ntLog["[prof] NumTrace TOTAL (", Length[diags], " diagrams): ", AbsoluteTime[] - ntT0, " s"];

  NTKernel[<|
    "Diagrams"  -> diags,
    "Env"       -> env,
    "NEnv"      -> nenv,
    "Frame"     -> frame,
    "Args"      -> args,
    "Dressings" -> dress
  |>]
];

(* One diagram -> {pure-scalar coeff, axis-id map, tensor components}. Components keep
   their factors un-expanded (heads, Plus-vertices, Times-structures); the recursive
   et compiler in Codegen turns Plus -> et::add, Times -> contract_all. *)
analyseDiagram[diagram_] := Module[{factors, tensorF, ids},
  factors = rewriteDressedNums @ splitSelfTraces[If[Head[diagram] === Times, List @@ diagram, {diagram}]];
  (* LAST CHANCE to close the fundamental-flavour deltas, and the last point at which promoting the
     residue into the SU(N) engine still works. rewriteDressedNums (just above) is what lifts a
     flavour delta out of an eager dressed numerator's Plus, so this is the first point at which a
     straddling chain is a flat product — and the partition below is the last point at which a
     promoted head can still be handed an axis id. A no-op unless the blind rules left something
     genuinely unclosable, so flows whose flavour lines close stay byte-identical. *)
  factors = promoteFlavResidue[factors];
  tensorF = Select[factors, ! scalarQ[#] &];
  (* Partition labels by sector: spinor (Dirac) axes get a disjoint high id range (>=100)
     so the et engine — which contracts axes by matching id — never fuses a spinor axis
     with a Lorentz/colour axis that happened to be numbered the same. *)
  ids     = With[{labs = DeleteDuplicates @ Flatten[allLabels /@ tensorF],
                  spn  = DeleteDuplicates @ Flatten[allSpinorLabels /@ tensorF]},
              With[{nonsp = DeleteCases[labs, Alternatives @@ spn]},
                Join[AssociationThread[nonsp -> Range[0, Length[nonsp] - 1]],
                     AssociationThread[spn -> Range[100, 99 + Length[spn]]]]]];
  <|
    (* contractFlavour collapses any flavour-δ chain that the dressing collection factored out of an
       eager numerator into the coeff (it straddled the Plus before, so FromFunKit could not). A no-op
       when collection is off / no flavour δ survives. *)
    "Coeff"      -> If[TrueQ[$ntDressCollect], contractFlavour[Times @@ Select[factors, scalarQ]],
                       Times @@ Select[factors, scalarQ]],
    "Ids"        -> ids,
    (* "Constant" is consumed by Codegen's per-component dispatch, where it means not merely
       "momentum-free" but "a constant SU(N) component": it is handed to compileColG, which only
       understands group heads. A momentum-free DIRAC structure (a bare closed spinor δ-loop, tr[1] = 4,
       with no ntVec/projector to carry momentum) is momentum-free yet carries no colour meaning, so
       compileColG emitted it as raw Mathematica (`colFacG[ntDeltaDirac[d2,d3], <|a1 -> 0, ...|>]`)
       into the generator .cpp — an instant clang failure. Reproduce with <P_2,T_2> of AqbqDirect8.
       Listing the Dirac heads here makes such a component NON-constant, routing it to
       splitColourGroupsInv/compileDirac like any other Dirac structure. NOTE this is only correct
       together with the collapsed-loop tr(1)=4 restoration in compileDirac (see nEmptyLoops there):
       orderDiracFacs drops δ connectors, so the loop arrives token-free and the runtime's split_loops
       would otherwise discard it, giving a kernel 4x too small.
       The SAME reasoning covers the pure-LORENTZ heads. A momentum-free closed metric loop
       (g_{μν} g^{μν} = D, the "ZAAqbq metric leak" shape) or a fully-contracted ntEpsilon pair from a
       γ5 trace is likewise momentum-free but not colour, and was emitted as
       `colFacG[ntMetric[v1,v2], <|v1 -> 0, ...|>]` straight into the generator .cpp. Reproduced with
       ntMetric[v1,v2] ntMetric[v2,v1] × ntVec[q1,a1] ntVec[ql,a1]. compileTInv/builderInv already
       handle both heads, so making the component non-constant is all that is needed.
       In short: "Constant" must mean "a constant SU(N) component", so EVERY non-SU(N) tensor head
       belongs in this list — it is not merely a momentum test. *)
    "Components" -> (<|"Factors" -> orderFactors[#],
                       "Constant" -> FreeQ[#, _ntVec | _ntTransProj | _ntLongProj |
                                              _ntElectricProj | _ntMagneticProj | _ntDressedNum | _ntDiracSlot |
                                              _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac |
                                              _ntMetric | _ntEpsilon]|> &
                     /@ connectedComponents[tensorF])
  |>
];

(* ---- per-diagram label validation ------------------------------------------------
   Run on ONE diagram AFTER expandBridges (a flat Times whose factors are heads, eager
   single-sector Plus vertices, Powers and scalars). The invariant, per diagram:
     * a label occurring ONCE  is a free (external) index of the whole trace;
     * a label occurring TWICE is contracted — the et engine pairs the two axes by id;
     * a label occurring 3+ times is MALFORMED. analyseDiagram assigns ONE id per DISTINCT
       label, make_eplan (axplan.hpp) pairs axes first-match-wins between exactly TWO
       operands, and contract_all is a LEFT FOLD — so four axes sharing an id silently
       contract as two independent pairs. No crash, just a wrong number. That is why this
       must Abort rather than Message.
   A Plus factor needs care and is why the old whole-net Cases[...,Infinity] count was
   useless: an eager vertex sum's summands LEGITIMATELY repeat the same free labels. A Plus
   is counted ONCE, via the free set its summands must agree on — that agreement is itself
   the et::add alignment precondition freeIdx[_Plus] assumes without checking. A summand's
   internal dummies are private and must not appear anywhere else. *)

NumTrace::flavleak = "a fundamental-flavour Kronecker delta survived BOTH contractFlavour passes \
and the promotion into the SU(N) engine. It is now neither contracted nor a tensor: scalarQ is a \
FreeQ over the nt* heads, so it counts as a scalar COEFFICIENT, its indices are invisible to \
labelsOf/freeIdx/checkLabels, and CForm would print it into the kernel as \
NumTracer_Private_flavDelta(F1, F2) — an undeclared identifier that GCC and Clang both parse, so \
the failure would surface as a link/compile error far from here, or compile silently wrong. Most \
likely cause: a flavour delta buried inside an eager (un-distributed) Plus that rewriteDressedNums \
did not lift out, so promoteFlavResidue saw it as already closed. Offending delta(s):\n`1`";

NumTrace::badlabel = "Diagram `1`: index label `2` occurs `3` times (expected 1 = free, \
2 = contracted). The et engine contracts axes by matching id, so `3` axes sharing this label \
are silently mis-paired into a wrong number. Offending diagram:\n`4`";
NumTrace::plusfree = "Diagram `1`: the summands of an eager (un-distributed) sum expose \
DIFFERENT free indices `2` — et::add cannot align them. Offending sum:\n`3`";
NumTrace::privclash = "Diagram `1`: label(s) `2` are private dummies of one factor but also \
occur outside it — a dummy-name collision between two independently generated objects. \
Offending diagram:\n`3`";

(* {exposed-multiset, private-set, bad-list} of a (sub)expression. *)
labelCensus[e_] := Which[
  tensorQ[e],
    With[{ls = labelsOf[e]},
      (* a label repeated WITHIN one head is a legal self-trace; splitSelfTraces resolves it *)
      {DeleteDuplicates[ls], Cases[Tally[ls], {l_, c_} /; c >= 2 :> l], {}}],

  Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 1 && ! scalarQ[e],
    (* n copies sharing the SAME labels: a closed self-contraction (see compileTInv) *)
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

  Head[e] === Plus,
    Module[{sub = labelCensus /@ (List @@ e), frees},
      frees = Sort /@ sub[[All, 1]];
      {If[frees === {}, {}, First[frees]],
       Union @@ sub[[All, 2]],
       Join[Join @@ sub[[All, 3]],
            If[Length[DeleteDuplicates[frees]] > 1, {{frees, "plus-free-mismatch"}}, {}]]}],

  True, {{}, {}, {}}];

(* Escape hatch: NT_NO_LABEL_CHECK=1 disables (the census is O(net), not a hot path). Anything
   falsy — unset, "", "0", "false" — leaves the check ON, which is the safe direction: this guard
   catches a label occurring more than twice, which otherwise becomes a silently wrong contraction. *)
$ntCheckLabels = !ntEnvFlag["NT_NO_LABEL_CHECK"];

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
