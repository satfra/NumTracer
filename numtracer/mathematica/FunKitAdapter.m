(* ::Package:: *)

(* FunKit adapter: rewrite a traced flow expression (FunKit's FormTracer head
   vocabulary) into the NumTracer DSL. A pure, local rewrite over a small closed
   set of heads — it does NOT touch the TensorBases internals.

   Expected input is a flow AFTER `// dressingRules // PropParam` (or SPParam): the
   dressings are ZA/Zc/RB expressions and the scalar dot products are already
   reduced to the runtime scalars (l1, p, cos1, ...). So only the TENSOR heads are
   rewritten here; their momenta are resolved by the frame in MakeKernel. (`sps`,
   the finite-T spatial scalar product, maps to ntSPS.)

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
     temporal rows. *)
  "transProjElectric" -> ntElectricProj, "transProjMagnetic" -> ntMagneticProj,
  "deltaLorentz" -> ntMetric,
  (* Dirac (spinor) sector. FunKit emits a slashed momentum as gamma[mu,..] * vec[q,mu]
     (the gamma carries the Lorentz axis), so ntGamma -> Dirac::gamma_axis. *)
  "gamma" -> ntGamma, "gamma5" -> ntGamma5, "deltaDirac" -> ntDeltaDirac,
  (* Charge conjugation. FormTracer's Dirac vocabulary is closed, so FunKit/TensorBases never
     trace C. A model supplies it as an INERT head on a hand-written vertex rule (or on a basis
     built with "Reduce"->False and "BuildProjectors"->False, which reaches no FORM call), and
     NumTracer gives it algebra. *)
  "ChargeConj" -> ntC,
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
  (* finite-T spatial scalar product, resolved by the frame as the SPATIAL dot (components 1..3).
     FunKit's component access vec[q, 0] (literal integer index) needs no rule here: it rides the
     "vec" -> ntVec map and DSL.m's integer-index classification routes it to q_0. *)
  "sps" -> ntSPS,
  (* finite-T spatial VECTOR vecs[q, mu] = {0, q_1, q_2, q_3}. Unlike sps this is a tensor LEG, so
     it gets a momentum of its own, ntSpatialVec[q] (DSL.m expandSpatialVecs / spatialVecFrame).
     Downstream a spatial slash vecs[q,mu] gamma[mu,d1,d2] is then an ordinary dslash. *)
  "vecs" -> (ntVec[ntSpatialVec[#1], #2] &)
|>;

(* The SU(N) group tokens, mapped to the N-parameterized heads with the group rank `n` injected
   as the leading argument. `nc` is the colour rank Global`Nc; `nf` the isospin rank (FromFunKit's
   "FlavourGroup" option). Colour structure (FCol/deltaAdjCol/TCol/deltaFundCol) builds against
   nc; the hand-rolled QM-model isospin tokens (the τ Yukawa generator, the pion f^{abc}
   self-coupling, the adjoint/fundamental isospin deltas) build against nf. Both go through the
   SAME four heads — the engine separates the groups by their disjoint contraction ids. *)
(* A generator with a FIXED (numeric) adjoint index, e.g. the Cartan direction TCol[3, i, j] of a
   Polyakov / A_0 background. SUNFac has no pinned adjoint index, and passing the literal through
   would make it an ordinary contraction LABEL that gets summed (tr(T^3 T^3) -> 4 instead of 1/2,
   silently). So pin it with a diagAdj that keeps only that one adjoint component.

   The index must be REMAPPED: FormTracer uses Gell-Mann ordering (SU(3) Cartans a = 3, 8), while
   NumTracer's generalized Gell-Mann generators put the N-1 diagonal ones LAST. Gell-Mann diagonal
   a = n^2-1 (n = 2..N) is 1-based NumTracer component N^2-N+n-1 (SU(3): 3 -> 7, 8 -> 8; SU(2):
   3 -> 3). A wrong remap is invisible to symmetric tests (tr(T^a T^a) = 1/2 for every a). Off-
   diagonal fixed indices are convention-dependent, so they are refused. *)
FromFunKit::cartan = "a generator with FIXED adjoint index `1` for SU(`2`): only the Cartan (diagonal) directions a = `3` can be pinned; the off-diagonal generators' positions in NumTracer's generalized Gell-Mann ordering are convention-dependent and would silently select the wrong one.";

ntCartanComponent[n_, a_] := Module[{m = Position[Table[k^2 - 1, {k, 2, n}], a]},
  If[m === {},
    Message[FromFunKit::cartan, a, n, Table[k^2 - 1, {k, 2, n}]];
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

(* FunKit's DECLARED head vocabulary (ShowFormTracerDefinitions[]): the closed-world list the guard
   in FromFunKit checks against. Every declared token must appear here, MAPPED OR NOT: a token
   absent from this list is invisible to the guard and leaks into the C++ as an opaque scalar. *)
$funKitHeads = {"FEx", "FTerm", "deltaLorentz", "vec", "vecs", "sp", "sps",
  "deltaDirac", "gamma", "gamma5", "ChargeConj", "sigma", "transProj", "longProj",
  "transProjElectric", "transProjMagnetic",
  "deltaAdjCol", "deltaFundCol", "FCol", "TCol", "epsAdjCol", "epsFundCol",
  "deltaAdjFlav", "deltaFundFlav", "fFlav", "tauFlav", "TFlav",
  "epsAdjFlav", "epsFundFlav", "deltaFlavFundGen", "epsLorentz"};
(* Handled outside `map`: TFlav by the hasIso rewrite in FromFunKit. *)
$ffHandledElsewhere = {"TFlav"};
FromFunKit::untranslated = "the FunKit head(s) `1` appear in the input but have no entry in \
$ffMap / sunMap. An untranslated head does NOT fail loudly downstream: DSL.m's scalarQ is a FreeQ \
over the KNOWN nt* heads, so an unknown head is classified as a SCALAR COEFFICIENT — its indices \
become invisible to labelsOf/freeIdx, the diagram reports spurious free (open) legs, checkLabels \
accepts it (open legs are legal), and the raw head is CForm'd into the generated C++. This is how \
epsFundCol/epsFundFlav went undetected. Add a $ffMap/sunMap entry, or refuse the input explicitly.";
(* epsLorentz is REFUSED on purpose (it falls through to FromFunKit::untranslated): it is the 3D
   SPATIAL epsilon, whereas ntEpsilon is 4D (DSL.m labelsOf, CodegenNets.m lorentzNetStr), so
   mapping one onto the other is a silent dimension error. Supporting it needs a spatial-delta head;
   ntMetric would wrongly include the temporal component. *)

FromFunKit::flavcount = "the fundamental-flavour sector would be closed against TWO different \
flavour counts in the same expression: the SU(N) engine uses rank `1` (the \"FlavourGroup\" option, \
defaulting to Nf when that is a bound integer and to 2 otherwise), while the blind \
contractFlavour folds a closed flavour loop to Nf = `2`. A single diagram can use both — a \
chain that closes cheaply beside a delta web the engine has to finish — so the coefficient would \
silently mix the two conventions rather than fail. Call SetNf[n] so Nf is the integer you \
mean, or pass \"FlavourGroup\" -> Nf explicitly.";

(* ---- adjoint Levi-Civita: SU(2) ONLY ----
   Only at rank 2 does the adjoint epsilon coincide with the structure constant, eps^{abc} = +f^{abc}
   (from NumTracer's T^a = sigma^a/2 and f^{abc} = -2i tr([T^a,T^b] T^c), sun_net.hpp; pinned by a
   test linear in the coefficient). Rewriting to ntSUNf also accepts an UNPAIRED epsilon. At rank N
   it carries N^2-1 indices; ShowFormTracerDefinitions[]' generic 3-index signature makes a 3-index
   epsAdjCol at Nc=3 look legal, so the refusal names that mistake. *)
FromFunKit::epsadj ="Adjoint Levi-Civita at SU(`1`) with `2` indices. NumTracer supports the \
adjoint epsilon ONLY at rank 2, where eps^abc coincides exactly with the structure constant f^abc \
(T^a = sigma^a/2, f = -2i tr([T^a,T^b]T^c); see sun_net.hpp:225) and is rewritten to ntSUNf[2,...]. \
That identification is SU(2)-SPECIFIC and does NOT generalise. At SU(`1`) the adjoint epsilon \
carries `3` indices, not `2` — ShowFormTracerDefinitions[] shows epsAdjCol[a,b,c] only as a generic \
three-index illustration, and taking that arity literally at rank > 2 is the likely mistake here. \
Aborting rather than guessing: a rank/arity mismatch would silently contract against the wrong \
group.";
adjEps[n_, idx__] := If[n === 2 && Length[{idx}] === 3, ntSUNf[2, idx],
  Message[FromFunKit::epsadj, n, Length[{idx}], n^2 - 1]; Abort[]];

(* ---- the ONE-ARGUMENT slash shorthand ----
   gamma[..., vec[p], ...] / gamma[..., vecs[p], ...] (no Lorentz index) is normally expanded by
   FormTracer itself, but can survive (finiteTenabled gating, FunKit's TRACY back-translation). Left
   alone it would put a nested head into ntGamma's Lorentz slot, failing far away, so expand it here
   first. Dispatch on head NAMES, as in the main map. *)
expandSlashShorthand[e_] := e //. (g_Symbol)[a___, (v_Symbol)[p_], b___] /;
    SymbolName[g] === "gamma" && MemberQ[{"vec", "vecs"}, SymbolName[v]] :>
  With[{mu = Unique["ffslash$"]}, v[p, mu] g[a, mu, b]];

(* Contract the flavour Kronecker deltas: their indices are disjoint from every tensor
   sector, so a chain collapses (delta[x,y] delta[y,z] -> delta[x,z]) and a closed loop
   delta[x,x] -> Nf. The result is a scalar power of Nf that the per-diagram coefficient
   carries (cancelling the projector's 1/Nf for a flavour-trivial flow like Zq).

   A Kronecker delta is SYMMETRIC, so all index orientations are matched explicitly; Orderless
   would say this in one rule but makes the //. matcher try permutations on a diagram-sized Times.

   The rules are SOUND but INCOMPLETE by design: each rewrite is an exact identity, and whatever
   is left (e.g. a delta WEB, not a chain) is handed to the SU(N) engine by promoteFlavResidue. *)
contractFlavour[e_] := e //. {
  flavDelta[x_, y_] flavDelta[y_, z_] :> flavDelta[x, z],
  flavDelta[y_, x_] flavDelta[y_, z_] :> flavDelta[x, z],
  flavDelta[x_, y_] flavDelta[z_, y_] :> flavDelta[x, z],
  flavDelta[x_, x_] :> Global`Nf,
  (* delta_{xy}^n = delta_{xy} for EVERY n >= 1, so summed over both indices it is Nf. *)
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

NumTrace::flavrankstale = "the flavour rank `1` handed over by the most recent FromFunKit call differs \
from Nf = `2`. A net that carries a fundamental-flavour delta always has rank Nf (FromFunKit::flavcount \
enforces it), so this net was built by an EARLIER FromFunKit call and a later one replaced the rank. \
Call NumTrace on a net right after the FromFunKit that built it.";

(* Hand the SU(N) engine whatever contractFlavour could not close.
   It runs from analyseDiagram, not FromFunKit, because contractFlavour runs TWICE: a chain that
   straddles an eager dressed numerator's Plus only closes after rewriteDressedNums lifts the common
   delta out. Promoting earlier would turn scalar Nf powers into SU(N) nets on every dressed quark flow.
   The no-residue path returns `factors` UNTOUCHED: analyseDiagram numbers axis ids by factor ORDER,
   so re-splicing would change the emitted code of unaffected flows. *)
promoteFlavResidue[factors_List] := Module[{flav, rest, closed, resid},
  If[FreeQ[factors, flavDelta], Return[factors]];
  flav   = Select[factors, ! FreeQ[#, flavDelta] &];
  closed = contractFlavour[Times @@ flav];
  resid  = DeleteDuplicates @ Cases[closed, flavDelta[__], {0, Infinity}];
  If[resid === {}, Return[factors]];
  If[! (IntegerQ[$ntFlavRank] && $ntFlavRank >= 1),
    Message[NumTrace::flavrank, Short[resid, 6]]; Abort[]];
  (* $ntFlavRank is a hand-off from the LAST FromFunKit call, not from the one that built this net *)
  If[IntegerQ[Global`Nf] && $ntFlavRank =!= Global`Nf,
    Message[NumTrace::flavrankstale, $ntFlavRank, Global`Nf]; Abort[]];
  (* A residue nested inside an eager Plus is promoted in place; the enclosing factor then fails
     scalarQ and correctly joins the tensor factors as an SU(N) Plus-vertex (compileColourSum). *)
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
FromFunKit[expr_, OptionsPattern[]] := Block[{$ntPlusMemo = <||>}, Module[{nf, map, hasIso, isoRewritten, res},
  nf  = OptionValue["FlavourGroup"] /. Automatic :> If[IntegerQ[Global`Nf], Global`Nf, 2];
  map = Join[$ffMap, sunMap[Global`Nc, nf]];
  (* ONE FLAVOUR COUNT: the blind contractFlavour closes a flavour loop to Global`Nf, the engine
     to `nf`. Both can contribute to one diagram, so a mismatch would silently mix conventions. *)
  If[! FreeQ[expr, (h_Symbol)[___] /; SymbolName[h] === "deltaFundFlav"] && nf =!= Global`Nf,
    Message[FromFunKit::flavcount, nf, Global`Nf]; Abort[]];
  (* ISOSPIN GENERATORS (quark-meson flows). The notebook symbol `TFlav` is the SU(nf) fundamental
     generator, TFlav[a,f1,f2] = (T^a)_{f1 f2}, with the singlet TFlav[0,f1,f2] = delta/Sqrt[2 Nf].
     When present, route the WHOLE fundamental-flavour sector into the SU(nf) engine (TFlav and the
     connecting deltas), so tr(T^a ... T^a) contracts; otherwise TFlav would leak as an opaque
     scalar. Gated on TFlav, so flavour-blind flows are byte-identical. *)
  hasIso = ! FreeQ[expr, Global`TFlav];
  If[hasIso, map["deltaFundFlav"] = (ntSUNDeltaFund[nf, ##] &)];
  isoRewritten = expandSlashShorthand @ If[hasIso,
    expr //. {Global`TFlav[0, f1_, f2_] :> ntSUNDeltaFund[nf, f1, f2]/Sqrt[2 Global`Nf],
              Global`TFlav[a_, f1_, f2_]  :> ntSUNT[nf, a, f1, f2]},
    expr];
  (* UNTRANSLATED-HEAD GUARD: an unknown head is classified downstream as a scalar coefficient and
     CForm'd into the C++ without any error (see FromFunKit::untranslated), so refuse any declared
     FunKit head with no map entry. Read `map`, never $ffMap: the hasIso branch above amends it. *)
  With[{present = DeleteDuplicates @ Cases[isoRewritten, (h_Symbol)[___] :> SymbolName[h], {0, Infinity}]},
    With[{leftover = Complement[Intersection[present, $funKitHeads], Keys[map], $ffHandledElsewhere]},
      If[leftover =!= {}, Message[FromFunKit::untranslated, leftover]; Abort[]]]];
  $ntDressCollect = TrueQ[OptionValue["DressingCollection"]];
  $ntFromFunKitDressCollect = $ntDressCollect;   (* NumTrace warns if it is called with the other value *)
  $ntFlavRank     = nf;   (* consumed by promoteFlavResidue, from DSL.m's analyseDiagram *)
  (* Normalize fixed Lorentz components before expandBridges tests whether a finite-T spatial
     slash is a collectible dressed Dirac numerator. The work is bound to `res` outside ntLog's
     arguments (see ntExportCpp in CodegenCommon.m). *)
  With[{ntT = First@AbsoluteTiming[res = contractFlavour @ expandBridges @ expandFixedComponents[
      isoRewritten //. (h_Symbol)[a___] /; KeyExistsQ[map, SymbolName[h]] :> map[SymbolName[h]][a]]]},
    ntLog["[prof] FromFunKit (head rewrite + expandBridges): ", ntT, " s"]];
  res]];
