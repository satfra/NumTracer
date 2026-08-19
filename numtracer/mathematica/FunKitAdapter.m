(* ::Package:: *)

(* FunKit adapter: rewrite a traced flow expression (FunKit's FormTracer head
   vocabulary) into the NumTracer DSL. A pure, local rewrite over a small closed
   set of heads — it does NOT touch the TensorBases internals.

   Expected input is a flow AFTER `// dressingRules // PropParam` (or SPParam): the
   dressings are ZA/Zc/RB expressions and the scalar dot products are already
   reduced to the runtime scalars (l1, p, cos1, ...). So only the TENSOR heads are
   rewritten here; their momenta are resolved by the frame in MakeKernel. (`sps`,
   the finite-T scalar product, is its own symbol and left untouched — finite-T
   support comes later.)

       flow = traceExprcbc // dressingRules // PropParam;
       net  = FromFunKit[flow];

   Structure:  FEx[...] = sum of terms,  FTerm[...] = product of factors.

   We dispatch on the head's *name* (SymbolName), not the symbol itself: the FunKit
   heads live in TensorBases`/FunKit` contexts not necessarily on $ContextPath when
   this file loads, so a literal `transProj` pattern would be a different symbol. *)

(* The FunKit-token -> NumTracer-head map, EXCEPT the SU(N) group tokens: those carry no rank
   of their own, so FromFunKit injects the rank N (a bound integer) into the new N-parameterized
   heads — see $sunMap below. This part is rank-independent. *)
$ffMap = <|
  "FEx" -> Plus, "FTerm" -> Times, "NonCommutativeMultiply" -> Times,
  "transProj" -> ntTransProj, "longProj" -> ntLongProj, "vec" -> ntVec,
  (* finite-T transverse split. TensorBases and NumTracer agree on both the argument order
     (momentum, mu, nu) and the conventions: P_E = P_T - P_M = delta_{mu 0} delta_{nu 0}
     + qs_mu qs_nu/|q_vec|^2 - q_mu q_nu/q^2, P_M = delta_{ij} - q_i q_j/|q_vec|^2 with vanishing
     temporal rows. The engine side (labelsOf/momentumOf/needsInvQ/needsInvSQ in DSL.m,
     builderInv -> eproj/mproj in Codegen.m) has been complete and gated by compare_ftproj_num all
     along; only this table was missing, so every FunKit flow with an E/M-split gluon propagator
     hit the untranslated-head trap the ROOT-CLASS GUARD below describes. *)
  "transProjElectric" -> ntElectricProj, "transProjMagnetic" -> ntMagneticProj,
  "deltaLorentz" -> ntMetric,
  (* Dirac (spinor) sector. FunKit emits a slashed momentum as gamma[mu,..] * vec[q,mu]
     (the gamma carries the Lorentz axis), so ntGamma -> Dirac::gamma_axis. *)
  "gamma" -> ntGamma, "gamma5" -> ntGamma5, "deltaDirac" -> ntDeltaDirac,
  (* flavour-TRIVIAL Kronecker delta -> a private head, contracted to a power of Nf below (its
     dimension Nf is symbolic). Correct ONLY for a flavour-blind closed loop (delta^{ii} = Nf,
     e.g. Zq/ZA quark loops). A fundamental flavour delta sitting INSIDE a τ-trace must instead
     be emitted as deltaFlavFundGen (-> ntSUNT/ntSUNDeltaFund via $sunMap) so it stays in the
     engine and does not break the trace by collapsing. *)
  "deltaFundFlav" -> flavDelta,
  (* any scalar product the notebook's reduction leaves (e.g. external-external SP
     constants sp[p_i,p_j] that SPParam does not touch) -> ntSP, resolved by the frame
     as the Euclidean dot of the components. For Zc/ZA nothing survives PropParam, so
     this is a no-op there. *)
  "sp" -> ntSP,
  (* finite-T spatial scalar product. Maps to ntSPS, resolved by the frame as the SPATIAL dot
     (components 1..3) of the momenta. FunKit's component access vec[q, 0] (a literal integer index)
     rides the "vec" -> ntVec map above and is routed to the scalar temporal component q_0 by the
     integer-index classification in DSL.m — no special rule needed here. *)
  "sps" -> ntSPS,
  (* finite-T spatial VECTOR vecs[q, mu] = the spatial part of q as a 4-vector ({0, q_1, q_2, q_3}).
     Unlike sps — a scalar the frame folds into the coefficient — this is a genuine tensor LEG, so it
     needs a momentum of its own: ntSpatialVec[q], which NumTrace pushes through sums and hands to
     the frame with the temporal slot zeroed (DSL.m expandSpatialVecs / spatialVecFrame). The whole
     point is that nothing downstream has to know: a SPATIAL SLASH vecs[q,mu] gamma[mu,d1,d2] — the
     dominant use, and what FormTracer's gamma[..., vecs[q], ...] shorthand expands to
     (FormTracer.m:1362) — is then an ordinary dslash against a different momentum. *)
  "vecs" -> (ntVec[ntSpatialVec[#1], #2] &)
|>;

(* The SU(N) group tokens, mapped to the N-parameterized heads with the group rank `n` injected
   as the leading argument. `nc` is the colour rank Global`Nc; `nf` the isospin rank (FromFunKit's
   "FlavourGroup" option). Colour structure (FCol/deltaAdjCol/TCol/deltaFundCol) builds against
   nc; the hand-rolled QM-model isospin tokens (the τ Yukawa generator, the pion f^{abc}
   self-coupling, the adjoint/fundamental isospin deltas) build against nf. Both go through the
   SAME four heads — the engine separates the groups by their disjoint contraction ids. *)
(* A generator carrying a FIXED (numeric) adjoint index — the Cartan directions T^3 / T^8 of a
   Polyakov / A_0 background, written TCol[3, i, j] exactly as in the FormTracer models. SetNc[3]
   switches colour to FormTracer`SU3fundexplicit, so TensorBases contracts those correctly on the
   Wolfram side; but SUNFac has NO pinned adjoint index. Passing the literal straight through as
   ntSUNT[nc, 3, i, j] makes the 3 an ordinary contraction LABEL, and a repeated label is SUMMED:
   tr(T^3 T^3) would come out as Σ_a tr(T^a T^a) = 4 instead of 1/2 — silently, with no error.
   Pin it instead with a diagAdj that keeps only that one adjoint component.

   THE INDEX IS NOT THE SAME NUMBER IN BOTH ENGINES. TensorBases/FormTracer use the standard
   Gell-Mann ordering, where the two SU(3) Cartan generators are a = 3 and a = 8. NumTracer builds
   GENERALIZED Gell-Mann generators with the diagonal (Cartan) ones LAST — verified by printing the
   generator diagonals: 0-based gen[6] = diag(1/2,-1/2,0) is Gell-Mann T^3 and gen[7] =
   diag(1,1,-2)/(2 Sqrt[3]) is Gell-Mann T^8, while gen[0..5] are off-diagonal. So Gell-Mann a must
   be REMAPPED, or the pin silently selects an off-diagonal generator (a bug invisible to any test
   that sums symmetrically over components — tr(T^a T^a) = 1/2 for EVERY a).

   For SU(N) the diagonal generators are the last N-1, and the Gell-Mann-convention diagonal indices
   are a = n^2-1 for n = 2..N. Hence 1-based NumTracer component = N^2-N+n-1. SU(3): 3 -> 7, 8 -> 8.
   SU(2): 3 -> 3 (coincides). Any other fixed index is an off-diagonal generator whose position in
   the generalized ordering is convention-dependent, so refuse it rather than guess.
   Verified component-sensitively via tr(T^a D) against a diagonal D. *)
ntCartanComponent[n_, a_] := Module[{m = Position[Table[k^2 - 1, {k, 2, n}], a]},
  If[m === {},
    Print["[NumTracer] ERROR: FromFunKit got a generator with FIXED adjoint index ", a,
      " for SU(", n, "). Only the Cartan (diagonal) directions a = ", Table[k^2 - 1, {k, 2, n}],
      " can be pinned: the off-diagonal generators' positions in NumTracer's generalized",
      " Gell-Mann ordering are convention-dependent and would silently select the wrong one."];
    Abort[]];
  n^2 - n + (m[[1, 1]] + 1) - 1];

ntSUNTPinnable[n_][a_, i_, j_] :=
  If[IntegerQ[a],
    With[{lbl = Unique["ntCartan$"], comp = ntCartanComponent[n, a]},
      ntSUNT[n, lbl, i, j] ntSUNDiagAdj[n, lbl, lbl, {comp -> 1}]],
    ntSUNT[n, a, i, j]];

sunMap[nc_, nf_] := <|
  "FCol" -> (ntSUNf[nc, ##] &), "deltaAdjCol" -> (ntSUNDeltaAdj[nc, ##] &),
  "TCol" -> (ntSUNTPinnable[nc][##] &), "deltaFundCol" -> (ntSUNDeltaFund[nc, ##] &),
  "fFlav" -> (ntSUNf[nf, ##] &), "deltaAdjFlav" -> (ntSUNDeltaAdj[nf, ##] &),
  "tauFlav" -> (ntSUNTPinnable[nf][##] &), "deltaFlavFundGen" -> (ntSUNDeltaFund[nf, ##] &),
  (* FUNDAMENTAL Levi-Civita: N indices (colour SU(3) -> 3, isospin SU(2) -> 2). Contracted into
     Kronecker deltas by expandFundEps in DSL.m — no engine token. *)
  "epsFundCol" -> (ntEpsFund[nc, ##] &), "epsFundFlav" -> (ntEpsFund[nf, ##] &),
  (* ADJOINT Levi-Civita, SU(2) ONLY — see adjEps. *)
  "epsAdjCol" -> (adjEps[nc, ##] &), "epsAdjFlav" -> (adjEps[nf, ##] &)
|>;

(* ---- adjoint Levi-Civita: SU(2) ONLY ---------------------------------------------------------
   The adjoint of SU(N) has dimension N^2-1, so its epsilon carries N^2-1 indices: 3 for SU(2), but
   8 for SU(3). ShowFormTracerDefinitions[] displays `epsAdjCol[a, b, c]` with a generic THREE-index
   signature, which makes a 3-index adjoint epsilon at Nc=3 look legal. It is not — that is the
   likely mistake, so the refusal below names it.

   !! SU(2) ONLY !!  At rank 2 (and ONLY at rank 2) the adjoint epsilon coincides with the structure
   constant: eps^{abc} = f^{abc}, exactly, coefficient +1. This is NOT a general fact — for SU(3),
   f^{abc} is not an epsilon at all (it is not totally antisymmetric in the same sense and has
   different nonzero entries), so nothing here may be read as an SU(N) statement.

   The +1 is derived from NumTracer's OWN conventions, not assumed: sun_net.hpp:233-247 builds the
   N=2 generators as T^a = sigma^a/2 in (x,y,z) order, and sun_net.hpp:225 defines
   f^{abc} = -2i tr([T^a,T^b] T^c). With [T^a,T^b] = i eps^{abc} T^c and tr(T^c T^d) = delta^{cd}/2
   this gives f^{abd} = eps^{abd}. It is still pinned by a test (a contraction LINEAR in the
   coefficient, so a sign flip cannot hide) rather than trusted.

   Rewriting to ntSUNf reuses an already-tested primitive and — unlike a pair-contraction scheme —
   accepts an UNPAIRED adjoint epsilon (eps^{abc} T^b T^c and friends), which is common and
   perfectly well-defined. *)
(* FunKit's DECLARED head vocabulary (ShowFormTracerDefinitions[]). The closed-world list the guard
   in FromFunKit checks against: every one of these must be mapped, or refused on purpose. *)
(* NOTE `vecs` is listed even though it IS mapped. That is the point of the list: it is the
   closed-world guard, and a token that is absent from it is not merely unmapped but INVISIBLE to
   the guard below — which is exactly how vecs used to reach the emitted C++ as an opaque scalar.
   Every finite-T token FormTracer declares (FormTracer.m:62) must appear here, mapped or refused. *)
$funKitHeads = {"FEx", "FTerm", "deltaLorentz", "vec", "vecs", "sp", "sps",
  "deltaDirac", "gamma", "gamma5", "sigma", "transProj", "longProj",
  "transProjElectric", "transProjMagnetic",
  "deltaAdjCol", "deltaFundCol", "FCol", "TCol", "epsAdjCol", "epsFundCol",
  "deltaAdjFlav", "deltaFundFlav", "fFlav", "tauFlav", "TFlav",
  "epsAdjFlav", "epsFundFlav", "deltaFlavFundGen", "epsLorentz"};
(* Handled outside `map`: TFlav by the hasIso rewrite above. *)
$ffHandledElsewhere = {"TFlav"};
FromFunKit::untranslated = "the FunKit head(s) `1` appear in the input but have no entry in \
$ffMap / sunMap. An untranslated head does NOT fail loudly downstream: DSL.m's scalarQ is a FreeQ \
over the KNOWN nt* heads, so an unknown head is classified as a SCALAR COEFFICIENT — its indices \
become invisible to labelsOf/freeIdx, the diagram reports spurious free (open) legs, checkLabels \
accepts it (open legs are legal), and the raw head is CForm'd into the generated C++. This is how \
epsFundCol/epsFundFlav went undetected. Add a $ffMap/sunMap entry, or refuse the input explicitly.";
(* epsLorentz is REFUSED on purpose, not mapped: it is the 3D SPATIAL epsilon (O(3) after the
   heat-bath split), whereas ntEpsilon is 4D and hard-wired to four labels (DSL.m labelsOf,
   Codegen.m builderInv and the Length[lst] == 4 reconstruction). Mapping one onto the other is a
   silent dimension error. The pair route is cheap (2 terms at D=3) but produces SPATIAL deltas, and
   NumTracer has no spatial-delta head — ntMetric is the 4D Euclidean metric and would wrongly
   include the temporal component. So the blocker is a missing head, not cost: the fix is a spatial
   delta, not more epsilon machinery. It falls through to FromFunKit::untranslated. *)

FromFunKit::flavcount = "the fundamental-flavour sector would be closed against TWO different \
flavour counts in the same expression: the SU(N) engine uses rank `1` (the \"FlavourGroup\" option, \
defaulting to Global`Nf when that is a bound integer and to 2 otherwise), while the blind \
contractFlavour folds a closed flavour loop to Global`Nf = `2`. A single diagram can use both — a \
chain that closes cheaply beside a delta web the engine has to finish — so the coefficient would \
silently mix the two conventions rather than fail. Call SetNf[n] so Global`Nf is the integer you \
mean, or pass \"FlavourGroup\" -> Global`Nf explicitly.";

FromFunKit::epsadj = "Adjoint Levi-Civita at SU(`1`) with `2` indices. NumTracer supports the \
adjoint epsilon ONLY at rank 2, where eps^abc coincides exactly with the structure constant f^abc \
(T^a = sigma^a/2, f = -2i tr([T^a,T^b]T^c); see sun_net.hpp:225) and is rewritten to ntSUNf[2,...]. \
That identification is SU(2)-SPECIFIC and does NOT generalise. At SU(`1`) the adjoint epsilon \
carries `3` indices, not `2` — ShowFormTracerDefinitions[] shows epsAdjCol[a,b,c] only as a generic \
three-index illustration, and taking that arity literally at rank > 2 is the likely mistake here. \
Aborting rather than guessing: a rank/arity mismatch would silently contract against the wrong \
group.";
adjEps[n_, idx__] := If[n === 2 && Length[{idx}] === 3, ntSUNf[2, idx],
  Message[FromFunKit::epsadj, n, Length[{idx}], n^2 - 1]; Abort[]];

(* ---- the ONE-ARGUMENT slash shorthand --------------------------------------------------------
   FormTracer accepts a slashed momentum written inside the gamma string — gamma[..., vecs[p], ...]
   and gamma[..., vec[p], ...], with no Lorentz index — and expands it itself (FormTracer.m:1362,
   1375). So this normally never reaches us. But the expansion is gated on finiteTenabled for the
   vecs form, and FunKit's TRACY back-translation (modules/TRACY/Tools.m:95) reconstructs the tokens
   textually from FORM output, so a shorthand CAN survive. Left alone it would hand ntGamma a nested
   head in a Lorentz slot: labelsOf returns that whole subexpression as a "label", and the failure
   surfaces far away. Expand it here, in the token vocabulary, before anything else looks.
   Dispatch on the head NAME for the same reason the main map does — the FunKit heads live in
   contexts that need not be on $ContextPath when this file loads. *)
expandSlashShorthand[e_] := e //. (g_Symbol)[a___, (v_Symbol)[p_], b___] /;
    SymbolName[g] === "gamma" && MemberQ[{"vec", "vecs"}, SymbolName[v]] :>
  With[{mu = Unique["ffslash$"]}, v[p, mu] g[a, mu, b]];

(* Contract the flavour Kronecker deltas: their indices are disjoint from every tensor
   sector, so a chain collapses (delta[x,y] delta[y,z] -> delta[x,z]) and a closed loop
   delta[x,x] -> Nf. The result is a scalar power of Nf that the per-diagram coefficient
   carries (cancelling the projector's 1/Nf for a flavour-trivial flow like Zq).

   A Kronecker delta is SYMMETRIC, so all four index orientations have to be matched. The
   head-to-tail rule alone leaves delta[x,y] delta[x,z] and delta[y,x] delta[z,x] untouched,
   so a chain that FunKit happened to emit in the other order silently failed to close.
   SetAttributes[flavDelta, Orderless] would say this in one rule, but it makes the //.
   matcher try argument permutations at every attempt on a diagram-sized Times — the blowup
   DSL.m:466 warns about — so the transposes are spelled out instead: cheap and deterministic.

   These rules are SOUND but INCOMPLETE, and that is by design. Each rewrite is an exact
   Kronecker identity, so stopping early is never wrong, only unfinished; whatever is left
   over is handed to the SU(N) engine by promoteFlavResidue below. A delta WEB (the boson
   tadpole's delta_ac delta_bd + delta_ad delta_bc against a loop delta) is not a chain and
   no local rewrite can close it — that is the residue's reason to exist. *)
contractFlavour[e_] := e //. {
  flavDelta[x_, y_] flavDelta[y_, z_] :> flavDelta[x, z],
  flavDelta[y_, x_] flavDelta[y_, z_] :> flavDelta[x, z],
  flavDelta[x_, y_] flavDelta[z_, y_] :> flavDelta[x, z],
  flavDelta[x_, x_] :> Global`Nf,
  (* delta_{xy}^n = delta_{xy} for EVERY n >= 1, so summed over both indices it is Nf. The rule
     used to match n === 2 only, so a cube survived and then read as a 3x-repeated label. *)
  Power[flavDelta[x_, y_], n_Integer /; n >= 2] :> Global`Nf
};

(* The SU(nf) rank FromFunKit routed the fundamental-flavour sector against, published for
   promoteFlavResidue (which runs later, from DSL.m's analyseDiagram). Automatic = no FunKit
   input in this session, i.e. a hand-built DSL net, which cannot contain flavDelta at all. *)
$ntFlavRank = Automatic;

NumTrace::flavrank = "a fundamental-flavour Kronecker delta survived the blind contraction and \
must be handed to the SU(N) engine, but no flavour rank is available (FromFunKit was never run, \
so $ntFlavRank is unset). This should be unreachable: flavDelta is produced by exactly one \
$ffMap entry. Offending delta(s):\n`1`";

(* Hand the SU(N) engine whatever contractFlavour could not close.
   WHY HERE AND NOT IN FromFunKit: contractFlavour runs TWICE. A flavour chain that straddles an
   eager dressed numerator's Plus is invisible to the first pass (FromFunKit), and only becomes a
   flat product once rewriteDressedNums has lifted the common delta out — which is why DSL.m runs
   it a second time. Measured on a qbq-shaped diagram: 3 deltas survive FromFunKit, and all 3 close
   after rewriteDressedNums. Promoting at the end of FromFunKit would therefore convert a scalar
   Nf power into an SU(N) net on every dressed quark flow (Zq/ZAqbq*/aqbq147/lambda3d/...), which
   is exactly the churn this design exists to avoid.

   The no-residue path returns `factors` UNTOUCHED rather than the contracted product: analyseDiagram
   numbers its axis ids from the ORDER of the tensor factors, so re-splicing the list would renumber
   them and change the emitted code for flows that are otherwise unaffected. So this is a strict
   no-op unless there is something the blind rules genuinely could not close. *)
promoteFlavResidue[factors_List] := Module[{flav, rest, closed, resid},
  If[FreeQ[factors, flavDelta], Return[factors]];
  flav   = Select[factors, ! FreeQ[#, flavDelta] &];
  closed = contractFlavour[Times @@ flav];
  resid  = DeleteDuplicates @ Cases[closed, flavDelta[__], {0, Infinity}];
  If[resid === {}, Return[factors]];
  If[! (IntegerQ[$ntFlavRank] && $ntFlavRank >= 1),
    Message[NumTrace::flavrank, Short[resid, 6]]; Abort[]];
  (* A residue nested inside an eager Plus is promoted in place; the enclosing factor then fails
     scalarQ and correctly joins the tensor factors as an SU(N) Plus-vertex (compileColGSum). *)
  closed = closed /. flavDelta[i_, j_] :> ntSUNDeltaFund[$ntFlavRank, i, j];
  rest   = Select[factors, FreeQ[#, flavDelta] &];
  Join[rest, If[Head[closed] === Times, List @@ closed, {closed}]]];

(* "FlavourGroup" -> Automatic resolves the isospin SU(N) rank to Global`Nf when that is a defined
   integer, else defaults to 2; an explicit integer overrides it. The colour rank is Global`Nc. *)
Options[FromFunKit] = {"FlavourGroup" -> Automatic, "DressingCollection" -> True};

(* Rewrite heads (injecting each SU(N) group's rank N into the new heads), then DISTRIBUTE
   (expandBridges turns the propagator-numerator structure sums into separate flat-product
   diagrams), THEN contract the flavour deltas — only once each diagram is a flat product do a
   flavour chain's links sit in one Times so they can collapse to a power of Nf (before
   distribution they straddle a Plus and cannot).
   "DressingCollection" -> True sets the gate BEFORE expandBridges (FromFunKit runs it before
   NumTrace), so dressed Dirac numerators are kept eager here too; pass the SAME value to NumTrace. *)
FromFunKit[expr_, OptionsPattern[]] := Module[{nf, map, hasIso, isoRewritten, res},
  nf  = OptionValue["FlavourGroup"] /. Automatic :> If[IntegerQ[Global`Nf], Global`Nf, 2];
  map = Join[$ffMap, sunMap[Global`Nc, nf]];
  (* ISOSPIN GENERATORS (quark-meson flows). The notebook auxiliary `TFlav` is the SU(nf) FUNDAMENTAL
     flavour generator: TFlav[a, f1, f2] = (T^a)_{f1 f2} for an adjoint index a != 0, and the flavour
     SINGLET TFlav[0,f1,f2] = deltaFundFlav[f1,f2]/Sqrt[2 Nf]. It is NOT in $ffMap (it is a notebook
     symbol, not a FunKit token), so without this it leaks through as an opaque scalar and the isospin
     trace never closes — leaving the external flavour indices dangling (the pion (-I) factors then fail
     to cancel and a spurious imaginary part survives). Route the WHOLE fundamental-flavour sector into
     the SU(nf) engine: map TFlav to ntSUNT / ntSUNDeltaFund, and switch the connecting fundamental
     deltas from the blind `flavDelta` (which only collapses a genuinely CLOSED line to Nf) to the
     in-engine ntSUNDeltaFund so the generator trace tr(T^a ... T^a) actually contracts. Gated on the
     presence of TFlav, so flavour-blind flows (Zq/ZA/ZAqbq1/4/7, ...) are byte-identical. *)
  (* ---- ONE FLAVOUR COUNT, NOT TWO -----------------------------------------------------------
     The two routes close a flavour loop against DIFFERENT numbers: the blind contractFlavour
     folds it to Global`Nf, the engine folds it to `nf` (the "FlavourGroup" option, which falls
     back to 2 when Global`Nf is not a bound integer). Both routes can contribute to ONE diagram
     — a chain that closes cheaply next to a web that the engine has to finish — so if the two
     disagree the coefficient silently mixes conventions. That is a wrong number, not a crash,
     and it is most likely exactly where it hurts: an Nf = 2+1 setup whose light group is SU(2).
     Refuse instead of guessing. Only checked when a fundamental-flavour delta is actually
     present; dispatch on the head NAME, as everywhere else in this file. *)
  If[! FreeQ[expr, (h_Symbol)[___] /; SymbolName[h] === "deltaFundFlav"] && nf =!= Global`Nf,
    Message[FromFunKit::flavcount, nf, Global`Nf]; Abort[]];
  hasIso = ! FreeQ[expr, Global`TFlav];
  If[hasIso, map["deltaFundFlav"] = (ntSUNDeltaFund[nf, ##] &)];
  isoRewritten = expandSlashShorthand @ If[hasIso,
    expr //. {Global`TFlav[0, f1_, f2_] :> ntSUNDeltaFund[nf, f1, f2]/Sqrt[2 Global`Nf],
              Global`TFlav[a_, f1_, f2_]  :> ntSUNT[nf, a, f1, f2]},
    expr];
  (* ---- THE ROOT-CLASS GUARD ------------------------------------------------------------------
     An untranslated FunKit head does NOT fail loudly downstream. DSL.m's scalarQ is a FreeQ over the
     KNOWN nt* heads, so an unknown head is classified as a SCALAR COEFFICIENT: its indices become
     invisible to labelsOf/freeIdx, the diagram reports spurious free (open) legs, checkLabels
     accepts it (open legs are legal), and the raw Mathematica head is CForm'd into the generated
     C++. That is exactly how epsFundCol/epsFundFlav went undetected through a whole debugging
     session — the visible symptom was 8 dangling indices three layers away from the cause.
     So: refuse any head from FunKit's declared vocabulary that has no entry in the map.
     NOTE the map is completed CONDITIONALLY just above (the hasIso branch promotes deltaFundFlav),
     so this must read `map`, never $ffMap — reading the wrong one is itself a way to be misled. *)
  With[{present = DeleteDuplicates @ Cases[isoRewritten, (h_Symbol)[___] :> SymbolName[h], {0, Infinity}]},
    With[{leftover = Complement[Intersection[present, $funKitHeads], Keys[map], $ffHandledElsewhere]},
      If[leftover =!= {}, Message[FromFunKit::untranslated, leftover]; Abort[]]]];
  $ntDressCollect = TrueQ[OptionValue["DressingCollection"]];
  $ntFlavRank     = nf;   (* consumed by promoteFlavResidue, from DSL.m's analyseDiagram *)
  ntLog["[prof] FromFunKit (head rewrite + expandBridges): ",
   (* Normalize fixed Lorentz components before expandBridges tests whether a
      finite-T spatial slash is a collectible dressed Dirac numerator. *)
   First@AbsoluteTiming[res = contractFlavour @ expandBridges @ expandFixedComponents[
      isoRewritten //. (h_Symbol)[a___] /; KeyExistsQ[map, SymbolName[h]] :> map[SymbolName[h]][a]]], " s"];
  res];
