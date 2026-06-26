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
  "sps" -> ntSPS
|>;

(* The SU(N) group tokens, mapped to the N-parameterized heads with the group rank `n` injected
   as the leading argument. `nc` is the colour rank Global`Nc; `nf` the isospin rank (FromFunKit's
   "FlavourGroup" option). Colour structure (FCol/deltaAdjCol/TCol/deltaFundCol) builds against
   nc; the hand-rolled QM-model isospin tokens (the τ Yukawa generator, the pion f^{abc}
   self-coupling, the adjoint/fundamental isospin deltas) build against nf. Both go through the
   SAME four heads — the engine separates the groups by their disjoint contraction ids. *)
sunMap[nc_, nf_] := <|
  "FCol" -> (ntSUNf[nc, ##] &), "deltaAdjCol" -> (ntSUNDeltaAdj[nc, ##] &),
  "TCol" -> (ntSUNT[nc, ##] &), "deltaFundCol" -> (ntSUNDeltaFund[nc, ##] &),
  "fFlav" -> (ntSUNf[nf, ##] &), "deltaAdjFlav" -> (ntSUNDeltaAdj[nf, ##] &),
  "tauFlav" -> (ntSUNT[nf, ##] &), "deltaFlavFundGen" -> (ntSUNDeltaFund[nf, ##] &)
|>;

(* Contract the flavour Kronecker deltas: their indices are disjoint from every tensor
   sector, so a chain collapses (delta[x,y] delta[y,z] -> delta[x,z]) and a closed loop
   delta[x,x] -> Nf. The result is a scalar power of Nf that the per-diagram coefficient
   carries (cancelling the projector's 1/Nf for a flavour-trivial flow like Zq). *)
contractFlavour[e_] := e //. {
  flavDelta[x_, y_] flavDelta[y_, z_] :> flavDelta[x, z],
  flavDelta[x_, x_] :> Global`Nf,
  Power[flavDelta[x_, y_], 2] :> Global`Nf
};

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
FromFunKit[expr_, OptionsPattern[]] := Module[{nf, map, hasIso, isoRewritten},
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
  hasIso = ! FreeQ[expr, Global`TFlav];
  If[hasIso, map["deltaFundFlav"] = (ntSUNDeltaFund[nf, ##] &)];
  isoRewritten = If[hasIso,
    expr //. {Global`TFlav[0, f1_, f2_] :> ntSUNDeltaFund[nf, f1, f2]/Sqrt[2 Global`Nf],
              Global`TFlav[a_, f1_, f2_]  :> ntSUNT[nf, a, f1, f2]},
    expr];
  $ntDressCollect = TrueQ[OptionValue["DressingCollection"]];
  contractFlavour @ expandBridges[
    isoRewritten //. (h_Symbol)[a___] /; KeyExistsQ[map, SymbolName[h]] :> map[SymbolName[h]][a]]];
