(* CodegenNets.m — net builders: lower one diagram component (Lorentz tensors, SU(N) colour/flavour
   factors, Dirac chains, dressed numerators, collected Dirac slots) to the C++ net literals
   (NetVal builders, SUNNet, DiracNet, DSlotOpt) that the emitted generator contracts numerically.
   Loaded by NumTracer.m via ntLoadPart, in the NumTracer`Private` context. *)

(* ---- Lorentz-factor emission ------------------------------------------------------------------
   Each tensor head becomes a call to one of the generator's wrapper templates
   (tproj / lproj / eproj / mproj / lmetric / lvec / leps / sc / contract / add), which the emitted
   generator defines over `NetVal`. `lorentzNetStr` emits the CALL TEXT for one head;
   `lorentzElemStr` emits the same head as the `Elem` aggregate that a collected Dirac slot's
   factor nets take instead. *)

lorentzNetStr[ntMetric[mu_, nu_], ids_, env_, nonzeroCompMask_] :=
  "lmetric<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ">()";

lorentzNetStr[ntVec[q_, mu_], ids_, env_, nonzeroCompMask_] :=
  "lvec<" <> ToString[ids[mu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[nonzeroCompMask[q]] <> ">()";

lorentzNetStr[ntTransProj[q_, mu_, nu_], ids_, env_, nonzeroCompMask_] :=
  "tproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[nonzeroCompMask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ">()";

lorentzNetStr[ntLongProj[q_, mu_, nu_], ids_, env_, nonzeroCompMask_] :=
  "lproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[nonzeroCompMask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ">()";

lorentzNetStr[ntMagneticProj[q_, mu_, nu_], ids_, env_, nonzeroCompMask_] :=
  "mproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[nonzeroCompMask[q]] <> ", " <> ToString[env[q]["InvS"]] <> ">()";

lorentzNetStr[ntElectricProj[q_, mu_, nu_], ids_, env_, nonzeroCompMask_] :=
  "eproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[nonzeroCompMask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ", " <> ToString[env[q]["InvS"]] <> ">()";

lorentzNetStr[ntEpsilon[a_, b_, c_, d_], ids_, env_, nonzeroCompMask_] :=
  "leps<" <> ToString[ids[a]] <> ", " <> ToString[ids[b]] <> ", " <> ToString[ids[c]] <> ", " <> ToString[ids[d]] <> ">()";

(* A Lorentz factor as a single `network::Elem{...}` literal (for a collected Dirac slot's per-option
   `netFacs`). Mirrors lorentzNetStr's id/momentum/atom resolution but emits the Elem aggregate the
   numeric backend appends to the net, rather than a NetVal builder. Field order (network.hpp):
   {kind, a, b, vid, inv, vlc, c, d, invS}. A projector's momentum rides `vid = env Base` (elem_to_nelem
   reconstructs it as {{1.0, vid}}); a vector's rides `vlc`. *)
lorentzElemStr[ntMetric[mu_, nu_], ids_, env_] :=
  "Elem{Elem::Metric, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", -1, -1, {}}";
lorentzElemStr[ntVec[q_, mu_], ids_, env_] :=
  "Elem{Elem::Vector, " <> ToString[ids[mu]] <> ", -1, -1, -1, {{1.0, " <> ToString[env[q]["Base"]] <> "}}}";
lorentzElemStr[ntTransProj[q_, mu_, nu_], ids_, env_] :=
  "Elem{Elem::ProjT, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[env[q]["Inv"]] <> ", {}}";
lorentzElemStr[ntLongProj[q_, mu_, nu_], ids_, env_] :=
  "Elem{Elem::ProjL, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[env[q]["Inv"]] <> ", {}}";
lorentzElemStr[ntMagneticProj[q_, mu_, nu_], ids_, env_] :=
  "Elem{Elem::ProjM, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", -1, {}, 0, 0, " <> ToString[env[q]["InvS"]] <> "}";
lorentzElemStr[ntElectricProj[q_, mu_, nu_], ids_, env_] :=
  "Elem{Elem::ProjE, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[env[q]["Inv"]] <> ", {}, 0, 0, " <> ToString[env[q]["InvS"]] <> "}";
lorentzElemStr[ntEpsilon[a_, b_, c_, d_], ids_, env_] :=
  "Elem{Elem::Epsilon, " <> ToString[ids[a]] <> ", " <> ToString[ids[b]] <> ", -1, -1, {}, " <> ToString[ids[c]] <> ", " <> ToString[ids[d]] <> "}";

scaleStr[str_, 1] := str;

scaleStr[str_, 1.] := str;

scaleStr[str_, s_] := "sc<numtracer::Lit<numtracer::Cx{" <> cppNum[s] <> ", 0.0}>>(" <> str <> ")";

(* no tensor factor at all (a pure-scalar product): no net, like compileLorentz's scalar branch *)
wrapContract[{}] := "";

wrapContract[{one_}] := one;

wrapContract[many_] := "contract(" <> StringRiffle[many, ", "] <> ")";

(* ---- memo keys for the net builders ----------------------------------------------------------
   The caches are keyed on the builders' TRUE argument, not the argument tuple as written:
   - `env` and `nonzeroCompMask` are fixed per generation (caches are cleared per generation), so they
     enter once, via the per-generation stamp $ctCtx, rather than being hashed on every call.
   - `ids` is the WHOLE diagram's label map, but the emitted string only sees labels resolved through
     it. So the key is `e` with its labels already substituted (ntCanonIds); keying on raw {e, ids}
     makes identical structures from different diagrams miss.
   - Momenta are dropped from the substitution (KeyDrop on the env keys): they resolve through
     `env`/`nonzeroCompMask`, and mapping them onto integers could conflate two different nets. *)

(* per-generation stamp of {env, nonzeroCompMask, frame}; set in mkGenerateKernel *)
$ctCtx = 0;

(* The substitution rules are constant per diagram, so the per-diagram loop precomputes them into
   $ntCanonRules (a Dispatch). The `===` guard is a correctness guard: if a caller passes a different
   diagram's `ids`, the stale rules would make the key collide with a different structure and the memo
   would silently return the wrong net. On a mismatch we rebuild the rules instead; on the hot path the
   comparison is pointer-identical. *)

$ntCanonIdsSrc = None;
$ntCanonRules = None;

ntCanonIds[e_, ids_, env_] :=
  e /. If[ids === $ntCanonIdsSrc,
         $ntCanonRules,
         Normal[KeyDrop[ids, Keys[env]]]];

(* Lorentz expression -> {netString, scalar}. MEMOISED (recursion included): the same projector
   structures recur across every diagram and branch, so calls far outnumber distinct arguments.
   $ctCache is cleared per generation in mkGenerateKernel. *)

compileLorentz[e_, ids_, env_, nonzeroCompMask_] := ntProfTimed["compileLorentz",
  With[{h = Hash[{ntCanonIds[e, ids, env], $ctCtx}]},
    Lookup[$ctCache, h, $ctCache[h] = compileLorentzBody[e, ids, env, nonzeroCompMask]]]];

compileLorentzBody[e_, ids_, env_, nonzeroCompMask_] := Which[
    tensorQ[e],
      {lorentzNetStr[e, ids, env, nonzeroCompMask], 1},
    Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 1 && !scalarQ[e],
    (* A TENSOR raised to an integer power = that many copies sharing the SAME index labels, i.e. a
       closed self-contraction (e.g. ntMetric[v1,v2]^2 = g_{v1 v2} g_{v1 v2} = D). Expand into n
       contracted copies so the numeric index elimination folds it to a number; otherwise it would reach
       the leak guard below. *)
      Module[{cs = Table[compileLorentz[e[[1]], ids, env, nonzeroCompMask], {e[[2]]}]},
        {wrapContract[cs[[All, 1]]], Times @@ cs[[All, 2]]}],
    Head[e] === Times,
      Module[{parts = List @@ e, sc, tn, cs},
        sc = Select[parts, scalarQ];
        tn = Select[parts, !scalarQ[#]&];
        cs = compileLorentz[#, ids, env, nonzeroCompMask]& /@ orderFactors[tn];
        {wrapContract[cs[[All, 1]]], (Times @@ sc) (Times @@ cs[[All, 2]])}],
    Head[e] === Plus,
      Module[{cs = compileLorentz[#, ids, env, nonzeroCompMask]& /@ (List @@ e)},
        If[!AllTrue[cs[[All, 2]], NumericQ],
          Message[MakeNTKernel::eagernn, e];
          Abort[]];
        {"add(" <> StringRiffle[MapThread[scaleStr, {cs[[All, 1]], cs[[All, 2]]}], ", "] <> ")", 1}],
    (* a genuine scalar coefficient: no builder, the expression IS the scalar *)
    scalarQ[e],
      {"", e},
    (* Anything else still carries a tensor head the branches above cannot lower; emitting it as a
       bare C++ scalar would silently drop its indices. Fail loudly instead. *)
    True,
      Message[MakeNTKernel::tleak, Short[e, 6]];
      Abort[]];

(* ---- colour factor patterns (used by splitColourGroups and compileColour) -------------------- *)

(* the SU(N) heads of $ntHeads that reach codegen *)
ctHeads = Blank @* ntHeadSym /@ Select[ntSectorHeads[{"Adjoint", "Fundamental"}], ! TrueQ[#["Rewritten"]] &];
$ctHeadPat = Alternatives @@ ctHeads;

colourEntangledQ[e_] := !FreeQ[e, $ctHeadPat];

(* A group head OR an integer power of one. The colour/rest split is a level-1 Cases/DeleteCases
   over the factor list, and a closed colour/flavour loop yields deltaFund[N,i,j]^2 (= N), a Power,
   not a head. Matched by the head alone it would stay in the Lorentz remainder and leak into the
   C++; collected here, compileColour expands it into repeated SUNNet factors. *)

ctFac = Alternatives[$ctHeadPat, Power[$ctHeadPat, _Integer?Positive]];

mergeColNet["SUNNet{}", b_] := b;

mergeColNet[a_, "SUNNet{}"] := a;

mergeColNet[a_, b_] :=
  "SUNNet{" <> StringDrop[StringDrop[a, 7], -1] <> "," <> StringDrop[StringDrop[b, 7], -1] <> "}";

(* ---- per-component diagonal-dressing registry (ntSUNDiag{Fund,Adj}) ------------------------
   A diag-dressed group δ dresses SELECTED components with distinctly-named scalar dressings (e.g.
   the Cartan directions of a condensate) and DROPS the rest. Its `spec` is a rules list
   {c1 -> expr1, …, Default -> defExpr}: `ci` are 1-based component indices, `expri` scalar dressing
   expressions; components with no rule (and no Default) vanish. colourFacStr registers each distinct
   dressing expression under a small integer id `dr` and bakes a per-component id vector
   (component → dr, -1 = drop) into the emitted sun<n>.diag{Fund,Adj}(...,{d0,…}) factor. The C++ seam
   folds the net to a SUNPoly over these ids, and the integrand multiplies in the runtime sum
   Σ_t coeff_t Π expr — an ordinary scalar-dressing token, no array. Reset per generation. *)

$diagDrTable = <||>; $diagDrByKey = <||>; $diagDrCounter = 0;

resetDiagDr[] := (
    $diagDrTable = <||>;
    $diagDrByKey = <||>;
    $diagDrCounter = 0;);

(* A component's dressing is a COMPLETE expression, kinematics included (e.g.
   `{1 -> Zu[scale], 2 -> Zd[scale]}`; each component may use its own scale). It is interned here and
   frame-resolved at emission; nothing is applied to it on the way out. *)

diagDrId[expr_] := Module[{key = expr},
    If[!KeyExistsQ[$diagDrByKey, key],
      $diagDrByKey[key] = $diagDrCounter;
      $diagDrTable[$diagDrCounter] = <|"Expr" -> expr|>;
      $diagDrCounter++];
    $diagDrByKey[key]];

(* Parse a diag-dressing spec into a per-component id vector of length `dim` (component 0..dim-1;
   1-based physics index = v+1). Named components get diagDrId[expr]; a Default -> expr rule
   fills the rest; unmatched components are -1 (dropped). *)

diagComp2Dr[spec_, dim_] := Module[{rules = Flatten[{spec}], named, def},
    named = Association[Cases[rules, (c_Integer -> nm_) :> (c -> diagDrId[nm])]];
    def = Cases[rules, (Default -> nm_) :> diagDrId[nm]];
    def =
      If[def === {},
        -1,
        First[def]];
    Table[Lookup[named, v + 1, def], {v, 0, dim - 1}]];

diagVecStr[vec_] := "{" <> StringRiffle[ToString /@ vec, ","] <> "}";

(* ---- scalar-dressing registry (symbolic dressing collection: ntDressedNum slots) -------------
   Each DISTINCT dressing atom — a maximal non-numeric multiplicative factor in a dressed numerator's
   per-structure coefficient (a propagator dressing `Zq[...]`, `hSigL[...]`, a regulator `RB[...]`, a
   composite denominator `Power[D,-1]`, a kinematic `cos1`) — is interned under a small integer id
   that keys the sub-term's dressing monomial and the generator's fm.dress table. So identical atoms
   across slots / diagrams share ONE `f[]` slot and the runtime evaluates each dressing call ONCE.
   The atom expression is already frame-resolved (ntSP/ntVec components substituted) so its C++ fill
   is `cppFlat[atom]`. Reset per generation alongside resetDiagDr. *)

$drTable = <||>; $drByKey = <||>; $drCounter = 0;

resetDr[] := (
    $drTable = <||>;
    $drByKey = <||>;
    $drCounter = 0;);

drAtomId[atom_] := Module[{key = atom},
    If[!KeyExistsQ[$drByKey, key],
      $drByKey[key] = $drCounter;
      $drTable[$drCounter] = atom;
      $drCounter++];
    $drByKey[key]];

(* Decompose a (frame-resolved) dressed-structure coefficient into {Cx numeric, {drAtomId…}}: numbers
   fold into the complex coefficient; a positive-integer power b^n expands to n atom copies; every other
   non-numeric factor is one atom (interned via drAtomId). *)

drDecompose[coeff_] := Module[{
    factors =
      If[Head[coeff] === Times,
        List @@ coeff,
        {coeff}],
    num = 1,
    ids = {}},
    Do[
      Which[
        NumberQ[f],
          num *= f,
        MatchQ[f, Power[_, _Integer?Positive]],
          ids = Join[ids, ConstantArray[drAtomId[First[f]], Last[f]]],
        True,
          AppendTo[ids, drAtomId[f]]],
      {f, factors}];
    {N[num, 17], Sort[ids]}];

(* Chunk a Lorentz polynomial into several inv nets of <= $ntInvChunk top-level terms each, so no
   single generated net-builder function becomes a giant nested add() that blows up the g++ -O0
   compile (a quark box with the full quark-gluon vertex basis can be tens of thousands of nodes in
   ONE net). Returns a list of {lorentzNetString, scalar}; the per-diagram combination sums the
   chunks (same colour + coeff) into one trace, so chunking is transparent. Zero -> dropped; a pure
   number -> a constant net (konst). *)

$ntInvChunk = 150;

chunkLorentz[lorExpr_, ids_, env_, nonzeroCompMask_] := Which[
    lorExpr === 0,
      {},
    scalarQ[lorExpr],
      {{"konst(" <> cppNum[lorExpr] <> ")", 1}},
    True,
      Module[{
        terms =
          If[Head[lorExpr] === Plus,
            List @@ lorExpr,
            {lorExpr}]},
        (compileLorentz[Total[#], ids, env, nonzeroCompMask]& /@ Partition[terms, UpTo[$ntInvChunk]])]];

(* ---- σ^{μν} folding: keep the bare γ-commutator as ONE token -----------------------------------
   The quark-gluon-vertex tensor σ^{μν}=(i/2)[γ^μ,γ^ν] arrives from FunKit as a bare 2-term
   antisymmetric γ-pair `Plus`  s·(γ(X)γ(Y) − γ(Y)γ(X))  (FunKit has no σ primitive). Left alone,
   splitColourGroups's `Expand` distributes it into TWO full Dirac traces. `foldDiracSigma`
   collapses that Plus into a single `ntSigma[legA, legB, din, dout]` token (each leg a slashed
   momentum {"slash",mom} or a free gluon id {"free",mu}), so the commutator is traced ONCE (the C++
   engine folds [A,B] as a block-diagonal 2×2 factor — `dcomm*` in network/dirac.hpp). The i/2 and any
   sign live in the SCALAR (the Plus is already a bare bracket). It is a pure OPTIMIZATION: when the
   Plus is not an UNAMBIGUOUS commutator the recognizer returns $Failed and it distributes. *)

(* the (type,value) leg of a gamma `g` within a term's factor list `tf`. The gamma's Lorentz label μ
   is carried by the factor(s) holding ntVec[_,μ): either a single ntVec[mom,μ] OR a nested momentum
   sum like (ntVec[l1,μ]−ntVec[p2,μ]) — the LOOP σ legs are momentum linear combinations (l1−p2, …),
   not single vecs. The leg is returned as a SORTED list of {coeff, q} pairs, each q the LITERAL ntVec
   momentum (an env key by construction: momentumOf collected every ntVec momentum into the env basis),
   so legStr maps it straight to a vlc without re-decomposing over a non-atomic basis. A μ carried by
   no ntVec ⇒ a free (open) gluon id {"free",μ}; an ill-formed leg (a μ-bearing factor that is not a
   plain coeff·ntVec sum) ⇒ $Failed (stay conservative — the commutator then distributes). *)

ntSigmaLeg[leg_, termFactors_] := Module[{mu = First[leg], momFactors, terms},
    momFactors = Times @@ Select[termFactors, !FreeQ[#, ntVec[_, mu]]&];
    If[momFactors === 1,
      Return[{"free", mu}]];
    terms =
      If[Head[momFactors] === Plus,
        List @@ momFactors,
        {momFactors}];
    Catch[
      {
        "slash",
        Sort[
          Function[trm,
              Module[{q, c},
                q = Cases[trm, ntVec[qq_, m_] /; m === mu :> qq, {0, Infinity}];
                c = trm /. ntVec[_, m_] /; m === mu :> 1;
                If[Length[q] =!= 1 || !FreeQ[c, ntVec] || !NumericQ[c],
                  Throw[$Failed, "sigleg"]];
                {c, First[q]}]
            ] /@ terms]},
      "sigleg"]];

(* chain-order the 2 gammas of one commutator term (spinor in->out) → {legFirst, legSecond, din, dout,
   scalar}, or $Failed if the term is not a clean isolated 2-γ pair. *)

ntSigmaTermInfo[term_] := Module[{tf, gs, g1, g2, first, second, din, dout},
    tf =
      If[Head[term] === Times,
        List @@ term,
        {term}];
    gs = Cases[tf, _ntGamma];
    If[Length[gs] =!= 2,
      Return[$Failed]];
    (* the two γ's (and their ntVecs) must be the term's ONLY tensor structure — no colour, projector,
       metric, other Dirac heads: anything else means the Plus is not a clean bare commutator. *)
    If[!FreeQ[tf, _ntGamma5 | _ntC | _ntDeltaDirac | _ntSigma | _ntSUNT | _ntSUNDeltaFund | _ntSUNf | _ntSUNDeltaAdj | _ntTransProj | _ntLongProj | _ntMetric | _ntEpsilon],
      Return[$Failed]];
    {g1, g2} = gs;
    Which[
      g1[[3]] === g2[[2]],
        first = g1;
        second = g2;
        din = g1[[2]];
        dout = g2[[3]],
      g2[[3]] === g1[[2]],
        first = g2;
        second = g1;
        din = g2[[2]];
        dout = g1[[3]],
      True,
        Return[$Failed]];
    {ntSigmaLeg[first, tf], ntSigmaLeg[second, tf], din, dout, Times @@ Select[tf, scalarQ]}];

(* {scalar, ntSigma[...]} for a 2-term Plus that is a bare γ-commutator, else $Failed. The two terms
   must share spinor endpoints, carry the SAME two legs in SWAPPED order, and have opposite scalar.
   Result uses term-1's leg order and scalar: c·([X,Y] block) (sign-symmetric — either term gives the
   same fold, since c_2 = −c_1 and swapping legs negates the bracket). *)

ntRecognizeComm[p_] := Module[{terms, a, b},
    If[Head[p] =!= Plus || Length[p] =!= 2,
      Return[$Failed]];
    terms = List @@ p;
    a = ntSigmaTermInfo[terms[[1]]];
    b = ntSigmaTermInfo[terms[[2]]];
    If[a === $Failed || b === $Failed,
      Return[$Failed]];
    (* same spinor endpoints *)
    If[a[[3]] =!= b[[3]] || a[[4]] =!= b[[4]],
      Return[$Failed]
    ];
    (* legs swapped between terms *)
    If[!(a[[1]] === b[[2]] && a[[2]] === b[[1]]),
      Return[$Failed]
    ];
    (* identical legs ⇒ [X,X]=0 *)
    If[a[[1]] === a[[2]],
      Return[$Failed]
    ];
    (* opposite scalar sign *)
    If[Simplify[a[[5]] + b[[5]]] =!= 0,
      Return[$Failed]
    ];
    {a[[5]], ntSigma[a[[1]], a[[2]], a[[3]], a[[4]]]}];

(* NT_NO_SIGMA_FOLD=1 disables the fold (the commutator then distributes into two traces): a safety
   switch and the baseline for measuring the fold's sub-term reduction. *)

$ntSigmaFold := !ntEnvFlag["NT_NO_SIGMA_FOLD"];

(* Fold every bare γ-commutator Plus among `factors` into a single ntSigma token (gated by
   $ntSigmaFold), recursively: find a factor that is a Plus recognisable as a commutator [A,B]
   (ntRecognizeComm), replace it with its two σ legs, recurse until none remain — so the antisymmetric
   γ-pair is never distributed into two separate Dirac traces. *)

foldDiracSigma[factors_List] := Module[{commPlus, recognized},
    If[!$ntSigmaFold,
      Return[factors]];
    commPlus = FirstCase[factors, q_ /; (Head[q] === Plus && ntRecognizeComm[q] =!= $Failed), Missing[], {1}];
    If[MissingQ[commPlus],
      Return[factors]];
    recognized = ntRecognizeComm[commPlus];
    foldDiracSigma[Join[DeleteCases[factors, commPlus, {1}, 1], {recognized[[1]], recognized[[2]]}]]];

(* ---- colour/Lorentz sector split for a colour-ENTANGLED Lorentz/Dirac component ---------------
   With the full quark-gluon vertex basis a quark-loop component is a SUM whose terms pair different
   colour orderings (T^{c1}T^{c2} vs T^{c2}T^{c1}) with different Dirac traces: colour and Dirac
   factor only per term, not globally. splitColourGroups therefore
     1. folds σ commutators, then EXPANDs only the entangled Pluses (those mixing colour with
        Dirac structure) into branches — the big pure-Lorentz angular polynomial stays factored;
     2. splits each branch into its colour product and a colour-free rest, compiled by compileDirac
        (gamma chain, traced numerically in C++) or chunkLorentz → {colourProduct, {{core, scal, restStr}…}};
     3. GatherBy colour product, emitting ONE {colourNet, bodyNets, scalar, restNets} entry per group,
        so the net count tracks the colour graph, not the branch×ordering explosion. The generator
        sums each group's branches at runtime. *)

(* a Plus that mixes colour with Dirac structure: expanded into branches *)
scgEntangledQ[x_] := Head[x] === Plus && (colourEntangledQ[x] || !FreeQ[x, _ntGamma | _ntGamma5 | _ntC | _ntDeltaDirac]);

splitColourGroups[factors0_, ids_, env_, nonzeroCompMask_] :=
  Module[{factors = foldDiracSigma[factors0], needExpand, keepAll, keepCol, keepRest, keepColLeakQ, keepDiracQ,
          distributed, terms, branchNets, groups},
    needExpand = Select[factors, scgEntangledQ];
    keepAll = Select[factors, !scgEntangledQ[#]&];
    (* keepAll is a factor of EVERY branch: split it into colour and remainder, and test the remainder,
       once per call rather than once per branch *)
    keepCol = Cases[keepAll, ctFac];
    keepRest = DeleteCases[keepAll, ctFac];
    keepColLeakQ = !FreeQ[keepRest, $ctHeadPat];
    keepDiracQ = !FreeQ[keepRest, $diracHeadPat];
    (* small: product of the entangled Pluses only *)
    distributed = Expand[Times @@ needExpand];
    terms =
      If[Head[distributed] === Plus,
        List @@ distributed,
        {distributed}];
    (* per branch -> {colourProduct, list-of-{bodyOrCore,scal,restStr}} *)
    branchNets =
      Function[term,
          Module[{termFactors = If[Head[term] === Times, List @@ term, {term}], termRest, colProd, rest},
            colProd = Times @@ Join[Cases[termFactors, ctFac], keepCol];
            termRest = DeleteCases[termFactors, ctFac];
            (* gammas + Lorentz + numeric coeff (no colour) *)
            rest = Join[termRest, keepRest];
            (* Level-1 DeleteCases only strips bare colour factors; one buried in an unexpanded Plus
               would leak into the C++ via lorentzNetStr. Fail here, where the offender is identifiable. *)
            If[keepColLeakQ || !FreeQ[termRest, $ctHeadPat],
              Message[MakeNTKernel::colrest, Short[DeleteDuplicates @ Cases[rest, $ctHeadPat, {0, Infinity}], 6], Short[rest, 8]];
              Abort[]];
            {
              colProd,
              If[keepDiracQ || !FreeQ[termRest, $diracHeadPat],
                {compileDirac[rest, ids, env, nonzeroCompMask]},
                (* gamma chain: {core, scal, projectorRest} *)
                ({#[[1]], #[[2]], ""}&) /@ chunkLorentz[Times @@ rest, ids, env, nonzeroCompMask]]}]
        ] /@ terms;
    groups = GatherBy[branchNets, First];
    (* A branch's `core` is a bare DiracNet contracted against its Lorentz rest only at runtime
       (numeric_value_netval), so branches cannot be fused into one net. A group carries them as
       parallel lists: bodyNet = {core_b…} (DiracNet, or a Lorentz NetVal for a gamma-free branch),
       restNet = {{rest_b, scal_b}…}; the generator sums Σ_b scal_b·numeric_value_netval(dnet_b, lnet_b). *)
    Function[group,
        Module[{colProd = group[[1, 1]], colNet, colScalar, branchRecs},
          {colNet, colScalar} =
            If[colProd === 1,
              {"SUNNet{}", 1},
              compileColour[colProd, ids]];
          (* each = {core, scal, restStr} *)
          branchRecs = Flatten[group[[All, 2]], 1];
          {colNet, branchRecs[[All, 1]], colScalar, ({#[[3]], #[[2]]}&) /@ branchRecs}]
      ] /@ groups];

(* ---- colour/group dialect -----------------------------------------------------------------------
   A constant SU(N) component is a product of structure constants, generators and Kronecker deltas,
   emitted as a `SUNNet` literal for the generator's numeric SU(N) contraction. Each head carries its
   group rank N as leading argument, so one net can mix several groups (colour SU(Nc) ⊗ flavour
   SU(Nf)); sun_net.hpp's sun_value_cx contracts each rank separately and multiplies. Each factor is
   minted by the per-rank `sun<n>` SUNEnv (declared in the generator main, the analogue of
   LorentzEnv), so the rank is written once, not on every factor. *)

colourFacStr[ntSUNf[n_, a_, b_, c_], ids_] :=
  "sun" <> ToString[n] <> ".f(" <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> "," <> ToString[ids[c]] <> ")";

colourFacStr[ntSUNDeltaAdj[n_, a_, b_], ids_] :=
  "sun" <> ToString[n] <> ".deltaAdj(" <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> ")";

colourFacStr[ntSUNT[n_, a_, i_, j_], ids_] :=
  "sun" <> ToString[n] <> ".T(" <> ToString[ids[a]] <> "," <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> ")";

colourFacStr[ntSUNDeltaFund[n_, i_, j_], ids_] :=
  "sun" <> ToString[n] <> ".deltaFund(" <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> ")";

(* per-component diagonal dressings: parse the spec into a per-component dressing-id vector
   (component → dr, -1 = drop; 1-based physics indices) and emit a diag factor carrying it. *)

colourFacStr[ntSUNDiagFund[n_, i_, j_, spec_], ids_] :=
  "sun" <> ToString[n] <> ".diagFund(" <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> "," <> diagVecStr[diagComp2Dr[spec, n]] <> ")";

colourFacStr[ntSUNDiagAdj[n_, a_, b_, spec_], ids_] :=
  "sun" <> ToString[n] <> ".diagAdj(" <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> "," <> diagVecStr[diagComp2Dr[spec, n^2 - 1]] <> ")";

(* CATCH-ALL, and it must stay LAST: the six rules above are the only lowerable colour factors.
   Without it a non-matching factor returns unevaluated and StringRiffle would ToString it into the
   generator .cpp, surfacing as a C++ syntax error far from its cause. Same role as
   MakeNTKernel::tleak in compileLorentz. *)

colourFacStr[e_, _] := (
    Message[MakeNTKernel::colleak, Short[e, 6]];
    Abort[]);

(* Constant colour product -> {colourNetString, factoredScalar}. NOT memoised: like compileDirac
   (see below), caching it was measured as a net loss. *)

compileColour[e_, ids_] := Module[
    {
      parts =
        If[Head[e] === Times,
          List @@ e,
          {e}],
      sc,
      tn},
    (* A colour/flavour SUM raised to a power cannot be expanded by repetition: the copies would share
       index labels (see MakeNTKernel::colpow). Refuse before the rewrite below can do it. *)
    Cases[
      parts,
      Power[b_Plus, k_Integer?Positive] /; !scalarQ[b] :>
        (
          Message[MakeNTKernel::colpow, k, Short[b, 6]];
          Abort[])];
    (* A colour/flavour factor raised to an integer power (e.g. deltaAdjFlav^2 from a CLOSED meson
       loop -> the flavour trace N^2-1) is expanded into repeated SUNNet factors so sun_value
       contracts the shared indices; colourFacStr handles a single head, not Power[head,k].
       Restricted to a BARE group head: repeating a Plus or Times duplicates labels instead of
       closing a self-trace. *)
    parts = parts /. Power[b_, k_Integer?Positive] /; MatchQ[b, $ctHeadPat] :> Sequence @@ ConstantArray[b, k];
    sc = Select[parts, scalarQ];
    tn = Select[parts, !scalarQ[#]&];
    {"SUNNet{" <> StringRiffle[colourFacStr[#, ids]& /@ tn, ", "] <> "}", Times @@ sc}];

(* ---- a colour/flavour component that is a SUM ------------------------------------------------
   `SUNNet` is a flat PRODUCT with no sum node or coefficient (mergeColNet splices literals by string
   surgery), so a sum cannot live in one net. It need not: colour folds to a scalar and the generator
   already sums several nets sharing a group. A summed component therefore lowers to a LIST of
   {net, scalar} branches, the shape chunkLorentz returns. Given the product of all the diagram's
   constant components, one Expand does the whole cross-product. Reached only from the
   constant-component path; splitColourGroups's colProd is flat by construction. *)

$ntColSumMaxBranches = 4096;

MakeNTKernel::colsum = "compileColourSum: a constant colour/flavour component expands to `1` summed branches (limit `2`). Each branch costs one emitted SUNNet and one net record, so this would blow up the generator. Raise $ntColSumMaxBranches if the flow genuinely needs it.";

compileColourSum[e_, ids_] := Module[{
    terms =
      With[{x = Expand[e]},
        If[Head[x] === Plus,
          List @@ x,
          {x}]]},
    If[Length[terms] > $ntColSumMaxBranches,
      Message[MakeNTKernel::colsum, Length[terms], $ntColSumMaxBranches];
      Abort[]];
    compileColour[#, ids]& /@ terms];

(* ---- spinor slots of the Dirac heads (used by the spinor-loop walk below) ---------------------- *)

(* the IN (row) spinor slot of a Dirac head; its OUT (column) slot is the other one *)
diracIn[h_] := First[spinorLabelsHead[h]];

(* SPINOR-SLOT SYMMETRY. The spinor delta is symmetric, δ_{ab} = δ_{ba}: it may be traversed either
   way and has no orientation to record (which is why the walk below can be undirected). Every other
   Dirac head (γ, γ5, σ, a collected numerator or vertex slot) is an ordinary matrix whose row and
   column slots are NOT interchangeable; traversing it against the arrow means the network wants its
   TRANSPOSE, which `orderDiracFacs` records with an `ntTransposed` wrapper. *)
diracSpinorSymmetricQ[_ntDeltaDirac] := True;
diracSpinorSymmetricQ[_]             := False;

(* ---- Dirac trace in the C++ generator (network/dirac.hpp `dirac_value`) ----------------------
   The gamma chain is emitted as a `DiracNet` literal and the generator traces the closed spinor loop
   NUMERICALLY (4×4 matrix products) against the Lorentz rest, as colour is folded by
   `SUNNet`/`compileColour`. *)

orderDiracFacs::open = "the spinor-loop walk consumed `1` of `2` token-bearing Dirac factors — a spinor loop did not close, and emitting it would silently drop γ structure (a collapsed trace). Loop factors:\n`3`";

(* The token-bearing factors of ONE closed spinor loop in trace order; spinor-δ connectors carry no
   token but are followed. UNDIRECTED cycle walk: each Dirac factor is an edge between its two spinor
   labels, every label has degree 2, so the walk is deterministic. Seed at `First[facs]` entering on
   its `diracIn`, then always leave by the OTHER endpoint. This also closes loops through a symmetric
   external spinor-δ (e.g. a scalar external projector); the trace is cyclic, so the start is free. *)
orderDiracFacs[facs_] :=
  Module[{nodeFacs = Association[], cur = 1, prevLabel, out = {}, seen = {}, labels, exitLabel, nexts, nTok, revQ},
    Do[
      Module[{ls = spinorLabelsHead[facs[[i]]]},
        (nodeFacs[#] = Append[Lookup[nodeFacs, #, {}], i])& /@ ls],
      {i, Length[facs]}];
    (* enter First on its diracIn so it exits on diracOut *)
    prevLabel = diracIn[facs[[1]]];
    While[
      !MemberQ[seen, cur],
      AppendTo[seen, cur];
      (* ORIENTATION. Following the cycle IS a matrix product: entering a factor on its `diracIn` (row)
         and leaving on its `diracOut` (column) uses it as declared; entering on the other slot means
         the network wants it TRANSPOSED:
           sum_{l0..l(n-1)} M1[l0,l1] M2[l1,l2] .. Mn[l(n-1),l0]  =  tr(N1 N2 .. Nn),   Nk = Mk or Mk^T.
         A reversed factor is only MARKED here; the engine transposes it. Do not replace this by a sign
         rule: the untransposed walk differs from the network by a sign segment parity does not predict. *)
      revQ = ! diracSpinorSymmetricQ[facs[[cur]]] && prevLabel =!= diracIn[facs[[cur]]];
      If[MatchQ[facs[[cur]], _ntGamma | _ntGamma5 | _ntC | _ntSigma | _ntDressedNum | _ntDiracSlot],
        AppendTo[out, If[revQ, ntTransposed[facs[[cur]]], facs[[cur]]]]];
      labels = spinorLabelsHead[facs[[cur]]];
      (* the OTHER endpoint *)
      exitLabel = First[DeleteCases[labels, prevLabel], Missing[]];
      If[MissingQ[exitLabel],
        Break[]];
      nexts = Select[DeleteCases[Lookup[nodeFacs, exitLabel, {}], cur], !MemberQ[seen, #]&];
      If[nexts === {},
        Break[]];
      prevLabel = exitLabel;
      cur = First[nexts]];
    (* The loop must consume EVERY token-bearing factor (δ connectors carry none); a walk that fell
       short would silently drop γ structure (e.g. a lone surviving γ5 → tr(γ5)=0). *)
    nTok = Count[facs, _ntGamma | _ntGamma5 | _ntC | _ntSigma | _ntDressedNum | _ntDiracSlot];
    If[Length[out] =!= nTok,
      Message[orderDiracFacs::open, Length[out], nTok, facs];
      Abort[]];
    (* Mixed-orientation loops (e.g. qq diquarks) are valid: reversed factors are transposed by the
       engine (DFac::transposed, or dtrslot(k) for a slot). The transposition is checked by section J
       of tests/test_numeric_contract.cpp against a brute-force index contraction. *)
    out];

(* A component may contain SEVERAL independent closed spinor loops (a quark loop + the
   projection-closed external line, tied together only by gluon propagators). orderDiracFacs walks ONE
   loop and would drop the rest; partition the Dirac factors into spinor-connected groups first
   (two factors share a loop iff they share a spinor index), then order each loop. Returns a LIST of
   ordered token lists (one per loop) — a single-loop component yields a one-element list. The C++
   contraction traces each loop separately and contracts their shared gluon legs (DFac::LoopSep).
   MEMOISED: the same chains recur across diagrams and colour branches, so distinct inputs are far
   fewer than calls. Pure function of `facs`; $odCache is cleared per generation in mkGenerateKernel. *)

$odCache = <||>;

orderDiracLoops[facs_] := With[{h = Hash[facs]},
    Lookup[$odCache, h, $odCache[h] = orderDiracLoopsBody[facs]]];

orderDiracLoopsBody[facs_] := If[Length[facs] <= 1,
    {orderDiracFacs[facs]},
    Module[
      {sp = spinorLabelsHead /@ facs, edges, g, comps},
      (* Two factors are adjacent iff they SHARE a spinor label: bucket factors by label and read the
         edges off the buckets (O(n), not an O(n^2) pair scan). Sort+DeleteDuplicates fix the edge order
         lexicographically, and so ConnectedComponents' component order (byte-identity), and keep each
         pair once (two factors may share BOTH labels). *)
      edges = Sort @ DeleteDuplicates @ Flatten[Subsets[#, {2}]& /@ Values @ GroupBy[Flatten[Table[{l, i}, {i, Length[facs]}, {l, sp[[i]]}], 1], First -> Last], 1];
      g = Graph[Range[Length[facs]], UndirectedEdge @@@ edges];
      comps = ConnectedComponents[g];
      (* a single spinor loop is walked on `facs` as given; multi-loop components are split, with
         indices sorted to keep the input order within each loop (byte-identity). *)
      If[Length[comps] <= 1,
        {orderDiracFacs[facs]},
        (orderDiracFacs[facs[[Sort[#]]]])& /@ comps]]];

(* frame resolver for dressed-numerator option coefficients (ntSP/ntVec[q,i] -> components). Set in
   mkGenerateKernel to the diagram's resolveScale; Identity when the dressed path is inactive. *)

$ntDressResolve = Identity;

(* ---- Dirac chain tokens -> C++ ------------------------------------------------------------------
   Shared by compileDirac (closed loops), diracSlotStrBody (an open chain inside a collected slot) and
   dressedSlotStrBody. Every slash momentum is a literal ntVec momentum, hence an env key; a missing
   one aborts rather than print a Missing[...] into the generator. *)
$diracHeadPat = ntSectorPat[{"Dirac"}];

compileDirac::slottransposed = "a transposed token reached the collected-slot emitter, which cannot represent one (orderOpenChain never marks a transpose). Token:\n`1`";

compileDirac::diracleak = "un-handled Dirac structure in a non-Dirac factor: a dressed propagator-numerator sum was NEITHER distributed NOR collected into ntDressedNum, so the numeric backend would silently drop or leak its gamma structure (a collapsed trace or untranslated C++). This is a front-end collection gap (collectibleDiracSumQ rejected a sum that distributeQ also skipped). Offending factor(s):\n`1`";

compileDirac::badtok = "compileDirac: `1` is not a Dirac chain token (expected ntGamma, ntGamma5, ntC, ntSigma, a slot, or ntTransposed of one). Aborting rather than emit it as a gamma.";

compileDirac::envmiss = "`1` momentum `2` is absent from env `3`.";

envBaseStr[q_, env_, what_] := (
  If[! KeyExistsQ[env, q],
    Message[compileDirac::envmiss, what, q, Keys[env]];
    Abort[]];
  ToString[env[q]["Base"]]);

(* a linear combination of momenta {{c1,q1},{c2,q2},…} -> the C++ vlc {{c1,Base1},{c2,Base2},…} *)
vlcCpp[pairs_List, env_, what_] :=
  "{" <> StringRiffle[("{" <> cppNum[#[[1]]] <> "," <> envBaseStr[#[[2]], env, what] <> "}") & /@ pairs, ", "] <> "}";

(* an ntSigma leg: a free open Lorentz leg (its axis id) or a slashed leg (a vlc) *)
sigmaLegCpp[{"slash", pairs_List}, ids_, env_] := vlcCpp[pairs, env, "σ slash leg"];
sigmaLegCpp[{"free", mu_}, ids_, env_] := ToString[ids[mu]];
$sigmaBuilder = <|{"free", "free"} -> "dcomm(", {"slash", "slash"} -> "dcomm_ss(",
                  {"free", "slash"} -> "dcomm_fs(", {"slash", "free"} -> "dcomm_sf("|>;

(* the bare DFac of one fixed token. A γ whose Lorentz leg is contracted with an ntVec[q,μ] factor
   (vecOf: μ -> q) is a slash, otherwise a free leg. ntTransposed marks a factor the spinor walk
   traversed against its declared direction: the engine multiplies its transpose. *)
fixedTokCpp[ntTransposed[g_], vecOf_, ids_, env_] := "dtr(" <> fixedTokCpp[g, vecOf, ids, env] <> ")";
fixedTokCpp[_ntGamma5, __] := "dg5()";
fixedTokCpp[_ntC, __] := "dc()";
fixedTokCpp[ntSigma[a_, b_, __], vecOf_, ids_, env_] :=
  $sigmaBuilder[{a[[1]], b[[1]]}] <> sigmaLegCpp[a, ids, env] <> ", " <> sigmaLegCpp[b, ids, env] <> ")";
fixedTokCpp[ntGamma[mu_, _, _], vecOf_, ids_, env_] :=
  If[KeyExistsQ[vecOf, mu],
    "dslash({{1.0," <> envBaseStr[vecOf[mu], env, "slash"] <> "}})",
    "dgamma(" <> ToString[ids[mu]] <> ")"];
fixedTokCpp[g_, __] := (Message[compileDirac::badtok, g]; Abort[]);

(* one token of a closed chain. A slot (dressed numerator or collected Dirac slot) is a DChainTok
   referencing the slot list, which is Sow'n under "slot" in chain order; a transposed slot is
   dtrslot(k) (reversed at the dress_enumerate splice). In a dressed chain every fixed DFac is
   wrapped as dtfix(DFac), so a transposed one reads dtfix(dtr(…)). *)
chainTokCpp[g : (_ntDressedNum | _ntDiracSlot), dressed_, vecOf_, ids_, env_, mask_] := (
  Sow[If[Head[g] === ntDressedNum, dressedSlotStr[g, env], diracSlotStr[g, ids, env, mask]], "slot"];
  "dtslot(" <> ToString[$ntSlotN++] <> ")");
chainTokCpp[ntTransposed[g : (_ntDressedNum | _ntDiracSlot)], rest__] :=
  StringReplace[chainTokCpp[g, rest], StartOfString ~~ "dtslot(" -> "dtrslot("];
chainTokCpp[g_, dressed_, vecOf_, ids_, env_, mask_] :=
  If[dressed, "dtfix(" <> fixedTokCpp[g, vecOf, ids, env] <> ")", fixedTokCpp[g, vecOf, ids, env]];
$ntSlotN = 0;

(* one ntDressedNum's options. Each option's scalar coefficient is frame-resolved then split into a
   complex number × dressing atoms (drDecompose); a "slash" structure's momenta become a vlc, an
   "ident" structure is the spinor identity (no token).
   MEMOISED: the same propagator numerator recurs in every diagram and colour branch, and resolving one
   is expensive. Pure given `env` and $ntDressResolve, both fixed per generation; $dsCache is cleared
   per generation in mkGenerateKernel alongside $ctCache. *)

$dsCache = <||>;

dressedSlotStr[gf : ntDressedNum[_, _, _], env_] := With[{h = Hash[{gf, $ctCtx}]},
    Lookup[$dsCache, h, $dsCache[h] = dressedSlotStrBody[gf, env]]];

(* Returns one {structStr, num, dr} triple per option: the dressing-free "DSlotOpt{…}" string
   (coeff 1, no dress atoms), the numeric coefficient (folded into the sub-term scalar) and the
   dress-atom ids (the sub-term's dressing monomial), so the trace table dedups on structure alone.
   The generator (emitNumericGenerator) expands the Cartesian product of the chain's slots' options
   into one single-option sub-term per combination. *)
dressedSlotStrBody[ntDressedNum[opts_, _, _], env_] := Function[opt,
          Module[{num, dr, vlcStr},
            {num, dr} = drDecompose[$ntDressResolve[opt[[1]]]];
            vlcStr = If[opt[[2, 1]] === "slash", vlcCpp[opt[[2, 2]], env, "dressed slash"], ""];
            (* ident → empty toks; slash → one dslash token. netFacs is empty: a propagator
               numerator has no open leg. *)
            {"DSlotOpt{Cx{1,0}, {}, {" <>
              If[opt[[2, 1]] === "slash", "dslash(" <> vlcStr <> ")", ""] <> "}, {}}", num, dr}]
        ] /@ opts;

(* ---- general collected Dirac slot → C++ DSlot literal (any open-leg count) ----------------------
   An ntDiracSlot option keeps its WHOLE structure (Dirac chain × Lorentz-net factors); here we split
   each option into DSlotOpt{coeff, {dress}, {toks}, {netFacs}} — the Dirac chain as a token list
   (dgamma/dslash/dcomm/dg5, open legs = ids of the free Lorentz tokens) and the Lorentz factors as
   network::Elem literals (the gluon propagator/metric that closes the open leg). Every label already
   has an id and every momentum an env slot (allLabels/momentumOf recurse into the slot via
   Cases[Infinity]); internal legs (the γ↔projector bridge) keep their own distinct ids, closed within
   the option, so no fresh-id allocation is needed. *)

(* Order one OPEN spinor chain din→dout, returning the token-bearing factors (γ/γ5/σ) in chain order;
   δ connectors are followed but carry no token. Walk from the din endpoint (degree 1) along spinor
   adjacency — like orderDiracFacs but seeded at a KNOWN endpoint (the chain is open, not a cycle). *)
orderOpenChain::open = "a collected Dirac slot's chain `1` -> ... could not be walked: `2` of its `3` token-bearing factor(s) were reached. The chain's spinor labels do not connect din to its tokens, so emitting them in input order would multiply them in the wrong order. Factors:\n`4`";

orderOpenChain[facs_, din_] :=
  Module[{nodeFacs = <||>, cur, prevLabel = din, out = {}, seen = {}, exitLabel, nexts, start,
          nTok = Count[facs, _ntGamma | _ntGamma5 | _ntC | _ntSigma]},
    Do[(nodeFacs[#] = Append[Lookup[nodeFacs, #, {}], i]) & /@ spinorLabelsHead[facs[[i]]], {i, Length[facs]}];
    start = Lookup[nodeFacs, din, {}];
    If[start === {},
      If[nTok == 0, Return[{}]];   (* an option with no Dirac token: an empty chain *)
      Message[orderOpenChain::open, din, 0, nTok, Short[facs, 6]]; Abort[]];
    cur = First[start];
    While[! MemberQ[seen, cur],
      AppendTo[seen, cur];
      If[MatchQ[facs[[cur]], _ntGamma | _ntGamma5 | _ntC | _ntSigma], AppendTo[out, facs[[cur]]]];
      exitLabel = First[DeleteCases[spinorLabelsHead[facs[[cur]]], prevLabel], Missing[]];
      If[MissingQ[exitLabel], Break[]];
      nexts = Select[DeleteCases[Lookup[nodeFacs, exitLabel, {}], cur], ! MemberQ[seen, #] &];
      If[nexts === {}, Break[]];
      prevLabel = exitLabel; cur = First[nexts]];
    If[Length[out] =!= nTok, Message[orderOpenChain::open, din, Length[out], nTok, Short[facs, 6]]; Abort[]];
    out];

(* Same memo keying as compileLorentz (see the memo-key note there): labels resolved through `ids`,
   the generation-fixed context (incl. nonzeroCompMask) in $ctCtx. *)
$dslCache = <||>;
diracSlotStr[gf : ntDiracSlot[_, _, _, _], ids_, env_, nonzeroCompMask_] := With[{h = Hash[{ntCanonIds[gf, ids, env], $ctCtx}]},
    Lookup[$dslCache, h, $dslCache[h] = diracSlotStrBody[gf, ids, env, nonzeroCompMask]]];

(* Returns one {structStr, num, dr} triple per option, like dressedSlotStrBody. *)
diracSlotStrBody[ntDiracSlot[opts_, din_, dout_, legs_], ids_, env_, nonzeroCompMask_] := Function[opt,
        Module[{num, dr, facs, vecOf, gammaLegs, diracFacs, lorFacs, toks, netFacs},
          {num, dr} = drDecompose[$ntDressResolve[opt[[1]]]];
          facs = If[Head[opt[[2]]] === Times, List @@ opt[[2]], {opt[[2]]}];
          vecOf = Association[Reverse[Cases[facs, ntVec[q_, m_] :> (m -> q)]]];
          gammaLegs = Cases[facs, ntGamma[gm_, _, _] :> gm];
          diracFacs = Select[facs, MatchQ[#, _ntGamma | _ntGamma5 | _ntC | _ntSigma | _ntDeltaDirac] &];
          (* Lorentz-net factors = non-Dirac tensors that are NOT a slash-vec (a slash-vec's μ is a γ leg) *)
          lorFacs = Select[facs, (tensorQ[#] && ! MatchQ[#, _ntGamma | _ntGamma5 | _ntC | _ntSigma | _ntDeltaDirac] &&
                          ! MatchQ[#, ntVec[_, m_ /; MemberQ[gammaLegs, m]]]) &];
          (* A slot's chain runs din->dout by construction (NumTrace::slotorient refuses an ambiguous
             one), so orderOpenChain never marks a transpose and the slot emitter cannot express one. *)
          toks = Function[gf2,
            If[MatchQ[gf2, _ntTransposed],
              Message[compileDirac::slottransposed, gf2]; Abort[]];
            fixedTokCpp[gf2, vecOf, ids, env]] /@ orderOpenChain[diracFacs, din];
          netFacs = lorentzElemStr[#, ids, env] & /@ lorFacs;
          (* dressing-free DSlotOpt; num and dr travel separately, as in dressedSlotStrBody *)
          {"DSlotOpt{Cx{1,0}, {}, {" <> StringRiffle[toks, ", "] <> "}, {" <>
            StringRiffle[netFacs, ", "] <> "}}", num, dr}
        ]] /@ opts;

(* A component's factor list -> {coreStr, scalar, restStr}: the Dirac chain(s) and the Lorentz rest,
   kept SEPARATE so splitColourGroups can share one rest across a colour group and the contraction
   stays deferred to the generator. A component without Dirac heads is just its Lorentz net.
     - the Dirac heads are ordered into independent spinor loops (orderDiracLoops), each a token
       string joined by loop separators;
     - a chain with a dressed numerator or a collected slot is returned as
       ntDressedCore[std::vector<DChainTok>{…}, slots] (slots: per-slot option lists, expanded into
       single-option sub-terms by emitNumericGenerator), else as DiracNet{…};
     - the remaining factors must be Dirac-free (a dressed numerator sum that was neither distributed
       nor collected would otherwise lose or leak its γ structure) and compile through compileLorentz.
   Not memoised: canonicalised calls repeat only ~2x, so the key costs more than it saves (measured
   as a net loss). *)
compileDirac[factors_, ids_, env_, nonzeroCompMask_] :=
  ntProfTimed["compileDirac", compileDiracBody[factors, ids, env, nonzeroCompMask]];

compileDiracBody[factors_, ids_, env_, nonzeroCompMask_] := Module[
    {diracFacs, vecOf, loops, dressed, loopStrs, slots, slashVecs, restFacs, restCompiled},
    diracFacs = Cases[factors, $diracHeadPat];
    If[diracFacs === {}, Return[Append[compileLorentz[Times @@ factors, ids, env, nonzeroCompMask], ""]]];
    (* μ -> q of an ntVec[q,μ]; Reverse so a duplicate μ keeps the FIRST match *)
    vecOf = Association[Reverse[Cases[factors, ntVec[q_, m_] :> (m -> q)]]];
    loops = orderDiracLoops[diracFacs];
    dressed = ! FreeQ[diracFacs, _ntDressedNum | _ntDiracSlot];
    {loopStrs, slots} = Block[{$ntSlotN = 0},
      Reap[StringRiffle[chainTokCpp[#, dressed, vecOf, ids, env, nonzeroCompMask] & /@ #, ", "] & /@ loops, "slot"]];
    slots = Flatten[slots, 1];
    (* the ntVec factors absorbed into a chain γ as its slash *)
    slashVecs = Cases[Join @@ loops /. ntTransposed[g_] :> g,
                      ntGamma[mu_, _, _] /; KeyExistsQ[vecOf, mu] :> ntVec[vecOf[mu], mu]];
    restFacs = DeleteCases[factors, Alternatives @@ Join[diracFacs, slashVecs]];
    If[! FreeQ[restFacs, $diracHeadPat],
      Message[compileDirac::diracleak, Select[restFacs, ! FreeQ[#, $diracHeadPat] &]];
      Abort[]];
    restCompiled = compileLorentz[Times @@ restFacs, ids, env, nonzeroCompMask];
    (* A loop of δ connectors only (a closed spinor δ-loop) has no token, so its segment is empty and
       the runtime's split_loops drops it: its tr(1) = 4 must be restored here. The bare chain deletes
       the empty segments and multiplies by 4 per loop. The dressed chain keeps every separator (the
       runtime counts collapsed loops from them) and only drops the empty segments' text. *)
    If[dressed,
      {ntDressedCore["std::vector<DChainTok>{" <> StringRiffle[DeleteCases[Riffle[loopStrs, "dtfix(dloopsep())"], ""], ", "] <> "}", slots],
       restCompiled[[2]], restCompiled[[1]]},
      {"DiracNet{" <> StringRiffle[DeleteCases[loopStrs, ""], ", dloopsep(), "] <> "}",
       restCompiled[[2]] * 4^Count[loopStrs, ""], restCompiled[[1]]}]];
