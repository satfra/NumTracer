(* ---- numeric (matrix-product) backend: component table + user symbols (task #22) -------------
   The numeric backend has NO sp-invariant basis: it evaluates scalar products numerically from each
   momentum's 4 COMPONENTS. So instead of an sp-invariant decomposition it needs
   only: the polynomial variables (the free user symbols), each fundamental momentum's 4 components as
   MPoly-builder C++, and the C++ fill formula for each symbol (a kernel arg, or a derived symbol like
   sin1 = sqrt(1-cos1^2)). Composite momenta resolve by component arithmetic via resolveComponents,
   exactly like the frame path — the user "Components" association is just a frame whose entries are
   polynomial. *)
(* Polynomialise a frame whose loop direction carries Sqrt[1-cos^2]: introduce a sin symbol per angle
   and record its definition. A user-supplied (already polynomial) "Components" assoc is a no-op. *)
(* Unit-loop spec (the COMPACT numeric parametrisation that matches the sp invariant count). The
   naive frame fallback bloats because the LOOP momentum is decomposed into angle PRODUCTS
   (l1·cos1, l1·sin1·cos2, … — degree-3), so every l·p_i expands into many angle monomials and the
   degrees compound through a projector chain. Instead write the loop as MAGNITUDE × UNIT-DIRECTION:
   comp_μ = l1 · Uμ where Uμ are degree-1 symbols (fill computes Uμ = dirμ from the kernel angles).
   Then l·p_i = l1·p·Σ Uμ u_i^μ is degree-1 in Uμ (as compact as sp(l,p_i)), and the bare-loop
   denominator Σ(l1·Uμ)² = l1²·ΣUμ² collapses to the MONOMIAL l1² under the unit constraint ΣUμ²=1
   (returned as a unit `group`), so its 1/l² atom cancels — exactly like inv's `rel`. Externals depend
   only on `p` (kept as the numeric p-vector). The loop is the momentum whose components carry `magSym`. *)

(* SPATIAL VECTORS (ntSpatialVec[q], from FormTracer's vecs) are NOT independent momenta: their
   components are the parent's with the temporal slot zeroed. Minting them a unit group of their own
   would be numerically correct — defs gives the duplicate ntU$ symbols the parent's very values —
   but it hides the identity from reduce_units, so cross terms like l·l̄ never collapse against the
   parent's Σ U² = 1 and the traces grow for nothing. Hold them out of the classification and derive
   them from the finished parent components instead. With no spatial vector in the frame every step
   below is a no-op, so existing kernels are untouched (Select, not Complement, everywhere — the key
   ORDER fixes the ntU$ numbering, and Complement would sort it). *)
spatialVecKeysOf[frame_] :=
  Select[Keys[frame], MatchQ[#, _ntSpatialVec] && KeyExistsQ[frame, First[#]] &];

unitLoopFrameSpec[frame_, pSym_, magSym_] := Module[{defs = <||>, groups = {}, n = 0, nf, svKeys},
    svKeys = spatialVecKeysOf[frame];
    nf =
      Association @
        KeyValueMap[
          Function[{q, comps},
            q ->
              If[SubsetQ[{pSym}, Variables[comps]],
                comps,
                (* external: numeric p-vector *)
                Module[
                  {grp = {}, rewrittenComps},
                  (* loop: comp_μ = magSym · Uμ *)
                  rewrittenComps =
                    Table[
                      Module[{dir = Coefficient[comps[[mu]], magSym], s},
                        If[dir === 0,
                          0,
                          s = Symbol["ntU$" <> ToString[n++]];
                          defs[s] = dir;
                          AppendTo[grp, s];
                          magSym s]],
                      {mu, 1, 4}];
                  AppendTo[groups, grp];
                  rewrittenComps]]],
          KeyDrop[frame, svKeys]];
    nf = Join[nf, Association[(# -> ReplacePart[nf[First[#]], 1 -> 0]) & /@ svKeys]];
    {nf, defs, groups}];

(* General frame -> {polyFrame, defs}: replace the trig sub-expressions in the frame components with
   fresh symbols so the components become polynomials in (magnitudes, cos's, sin's). Introduces a
   symbol for each loop polar factor Sqrt[1-cos^2] (a sin) and each bare-angle Cos/Sin (the 4-point
   frame's φ); `defs` maps each new symbol back to its closed form (so the kernel can fill it). The
   fallback path when the unit-loop parametrisation (unitLoopFrameSpec) does not apply. *)

polyFrameSpec[frame_] := Module[{defs = <||>, vals, rules = {}, polyFrame},
    vals = PowerExpand[Values[frame]];
    (* loop polar factor Sqrt[1-cos^2] -> a sin symbol *)
    Do[
      With[{s = Symbol["ntSin$" <> SymbolName[c]]},
        defs[s] = Sqrt[1 - c^2];
        AppendTo[rules, Sqrt[1 - c^2] -> s]],
      {c, DeleteDuplicates @ Cases[vals, Power[1 - c_^2, 1/2] :> c, Infinity]}];
    (* trig of a bare angle symbol (e.g. Cos[phi], Sin[phi] in the 4-point frame) -> cos/sin symbols *)
    Do[
      With[{s = Symbol["ntCos$" <> SymbolName[a]]},
        defs[s] = Cos[a];
        AppendTo[rules, Cos[a] -> s]],
      {a, DeleteDuplicates @ Cases[vals, Cos[a_Symbol] :> a, Infinity]}];
    Do[
      With[{s = Symbol["ntSinA$" <> SymbolName[a]]},
        defs[s] = Sin[a];
        AppendTo[rules, Sin[a] -> s]],
      {a, DeleteDuplicates @ Cases[vals, Sin[a_Symbol] :> a, Infinity]}];
    polyFrame = vals /. rules;
(* GENERAL radicals. The rewrites above only catch the shapes the symmetric-point frames happen to
   produce: Sqrt[1-c^2] for a bare symbol c, and trig of a bare angle. An arbitrary algebraic frame
   produces radicals of COMPOSITE arguments -- e.g. the (S0,S1,SPhi) three-point frame, whose
   magnitudes are Sqrt[S0^2 (1 - S1 Sin[SPhi])] and whose Q-k opening angle brings a
   1/Sqrt[1 - S1^2 Sin[SPhi]^2]. Those survive to numericComponents, which then rejects the frame
   outright ("fractional power of a symbol remains") -- so any such frame was simply unusable.

   A radical is just a scalar function of the runtime coordinates, exactly like the sin symbols
   above, so give each remaining one a fresh symbol and record its closed form in `defs`; the kernel
   fills it via cppFlat, which emits the sqrt inline. One symbol per distinct (base, exponent) pair
   so that negative half-powers become plain symbols too rather than 1/poly, keeping the components
   polynomial as numericComponents requires.

   ReplaceAll (not ReplaceRepeated) is what makes nesting safe: it rewrites the OUTERMOST radical
   and does not descend into the replacement, so a nested radical stays inside the recorded closed
   form and is emitted as part of that one C++ expression. *)
    Module[{n = 1, radSeen = <||>},
      polyFrame =
        polyFrame /. (pw : Power[b_, _Rational] /; ! NumericQ[b]) :>
          (If[! KeyExistsQ[radSeen, pw],
             With[{new = Symbol["ntRad$" <> ToString[n++]]},
               radSeen[pw] = new;
(* Record the CLOSED FORM in the raw runtime coordinates, not in the symbols minted above. By this
   point the radical is expressed in ntSin$/ntCos$/ntSinA$, and `vfill` emits a def verbatim via
   cppFlat -- so leaving them in produces C++ that names another fill slot instead of computing it:
       f[11] = sqrt(1. - ntSinA$SPhi * S1);   // 'ntSinA$SPhi' was not declared in this scope
   Those symbols are fill slots, not variables, and only the ones appearing in the COMPONENTS get a
   slot at all. Resolving `defs` here keeps each definition self-contained. *)
               defs[new] = pw //. Normal[defs]]];
           radSeen[pw])];
    polyFrame = Association @ Thread[Keys[frame] -> polyFrame];
    {polyFrame, defs}];

(* MIXED unit-loop spec: loop momenta as magnitude × unit-direction, externals through the general
   polynomialiser. The unit-loop fast path (unitLoopFrameSpec) requires every external to depend
   only on `pSym`, so a frame whose externals carry shape coordinates — the (S0,S1,SPhi) two-
   momentum frames of the with_mesons lambda1L3D class — used to fall ALL the way to polyFrameSpec:
   the loop was expanded into cos1/sin1·cos2/… component products, every trace polynomial carried
   those products at high degree, the bare-loop denominator Σ(l1·Uμ)² never collapsed to the
   monomial l1² (so divThroughMonomialAtoms could not cancel the 1/l² atoms), and the emitted
   traces blew up ~70x against the same derivative at fast-path kinematics.

   But the LOOP is frame-independent — its components are ∝ l1 in every vacuum frame — so it can
   keep the exact unit-direction treatment (opaque degree-1 ntU$ symbols + the ΣUμ²=1 group)
   regardless of what the externals look like; only the externals need polyFrameSpec's
   radical-minting. That is what this spec does, per momentum.

   Measured warning (lambda3d_small fixture, 2026-08-07): the OTHER seemingly-obvious repair —
   keeping the polyFrame cos/sin symbols and handing reduce_units the `{cos, ntSin$cos}` pair
   groups — makes emitted traces 7.3x BIGGER, not smaller: sin²→1−cos² substitution expands
   monomials binomially at high even powers, and the l1² collapse it buys back is worth far less.
   The unit-direction form wins because the loop enters traces at degree 1 per component, not
   because of the trig identity. Do not resurrect the pair-group variant.

   NT_NO_UNIT_GROUPS: escape hatch back to the pre-2026-08-07 behaviour (whole frame through
   polyFrameSpec). Exists for grading — the lambda3d_small control kernel is generated under it,
   so the lever is validated numerically against the pristine path — and as the rollback. *)

(* SPATIAL unit loop (the finite-T case): the temporal component (slot 1 ↔ C++ component 0) is an
   independent magSym-free coordinate (l0 / a Matsubara frequency), while the SPATIAL components
   are all ∝ magSym with at least one nonzero direction. O(3) is intact at T > 0, so the spatial
   direction is a polar-parametrized UNIT 3-vector — exactly the property the unit-group rewrite
   Σ U² = 1 needs, now over the spatial components only. The bare-loop denominator then reduces to
   the TWO-TERM l0² + l1² (divThroughPolyAtoms' case) instead of an angle-product polynomial.
   Requires the frame to parametrize the spatial part polar-wise (dirs unit-norm), the same
   assumption the full unit-loop branch makes about all four components. *)
(* a frame component that is magSym times a magSym-free coefficient *)
magPropQ[c_, magSym_] := Simplify[c - Coefficient[c, magSym] magSym] === 0;
(* all four components of an (already PowerExpand-ed) frame vector are magSym-proportional *)
fullLoopQ[cc_, magSym_] := AllTrue[Range[4], magPropQ[cc[[#]], magSym]&];

unitLoopSpatialQ[comps_, magSym_] := Module[{cc = PowerExpand[comps]},
  FreeQ[cc[[1]], magSym] &&
  AllTrue[Range[2, 4], magPropQ[cc[[#]], magSym]&] &&
  AnyTrue[Range[2, 4], (Coefficient[cc[[#]], magSym] =!= 0)&]];

unitLoopMixedOkQ[frame_, magSym_] :=
  !ntEnvFlag["NT_NO_UNIT_GROUPS"] &&
    (* at least one momentum is a magSym-proportional loop — either FULLY (all four components,
       the vacuum case) or SPATIALLY (finite T: an independent temporal l0 rides along) — and NO
       momentum mixes magSym with other coordinates inside a component *)
    With[{cc = PowerExpand[Values[frame]]},
      AnyTrue[cc, (fullLoopQ[#, magSym] || unitLoopSpatialQ[#, magSym])&] &&
      AllTrue[cc, (fullLoopQ[#, magSym] || unitLoopSpatialQ[#, magSym] || FreeQ[#, magSym])&]];

(* SHARED DIRECTIONS (2026-09-24). Loop tags that are the same vector up to the temporal slot -- a
   finite-T frame's l1 and its fermionic partner lf1 = l1 + (pi T, 0) from frameShiftedLoop, or any
   two keys with identical direction coefficients -- used to mint one ntU$ set and one unit group EACH.
   The traces then carried two independent copies of the same three direction numbers, so monomials
   in the copies could neither merge (cos1 * cos1' never became cos1^2) nor reduce against the SAME
   group's Sum U^2 = 1, and every polynomial downstream grew: ZA4 had MPoly vars {f0, l1, U2..U4,
   U5..U7, p, T} with U5..U7 == U2..U4 numerically. Keys are now matched on their direction list, so
   the second tag reuses the first one's symbols and group. `dirSyms` memoises by (kind, directions). *)
unitLoopMixedFrameSpec[frame_, magSym_] := Module[
    {svKeys, loopKeys, spatKeys, extFrame, pf, defs, groups = {}, n = 0, nf, nfS, dirSeen = <||>, dirSyms},
    (* direction list -> list of ntU$ symbols (0 where the direction vanishes); one unit group per
       DISTINCT direction list, minted on first sight and reused by every later key that matches *)
    dirSyms[kind_, dirs_List] :=
      With[{key = {kind, Simplify /@ dirs}},
        If[KeyExistsQ[dirSeen, key],
          dirSeen[key],
          Module[{grp = {}, syms},
            syms =
              Map[
                Function[dir,
                  If[dir === 0,
                    0,
                    With[{s = Symbol["ntU$" <> ToString[n++]]},
                      defs[s] = dir;
                      AppendTo[grp, s];
                      s]]],
                dirs];
            AppendTo[groups, grp];
            dirSeen[key] = syms]]];
    (* spatial vectors are derived from their parent at the end, never classified — see
       spatialVecKeysOf. Held out of loopKeys/spatKeys/extFrame so they mint nothing of their own. *)
    svKeys   = spatialVecKeysOf[frame];
    loopKeys = Select[Keys[frame], !MemberQ[svKeys, #] && fullLoopQ[PowerExpand[frame[#]], magSym]&];
    spatKeys = Select[Keys[frame],
      (!MemberQ[svKeys, #] && !MemberQ[loopKeys, #] && unitLoopSpatialQ[frame[#], magSym])&];
    (* externals — AND the spatial loops' temporal components — go through polyFrameSpec's
       sin/cos/radical minting, verbatim: a spatial loop is handed in as {l0, 0, 0, 0} so its
       temporal coordinate gets exactly the same treatment an external's would. *)
    extFrame = Join[
      KeyDrop[frame, Join[loopKeys, spatKeys, svKeys]],
      Association @ Map[# -> {PowerExpand[frame[#]][[1]], 0, 0, 0}&, spatKeys]];
    {pf, defs} = polyFrameSpec[extFrame];
    (* full loops: comp_μ = magSym · ntU$n, one unit group per loop — unitLoopFrameSpec's treatment *)
    nf =
      Association @
        Map[
          Function[q,
            q ->
              Module[{cc = PowerExpand[frame[q]], syms},
                syms = dirSyms["full", Table[Coefficient[cc[[mu]], magSym], {mu, 1, 4}]];
                magSym syms]],
          loopKeys];
    (* spatial loops (finite T): temporal component = polyFrameSpec's minted form, spatial
       components = magSym · ntU$n with the unit group over the SPATIAL directions only
       (Σ_{i=1..3} U_i² = 1 — O(3) polar parametrization). *)
    nfS =
      Association @
        Map[
          Function[q,
            q ->
              Module[{cc = PowerExpand[frame[q]], syms},
                syms = dirSyms["spatial", Table[Coefficient[cc[[mu]], magSym], {mu, 2, 4}]];
                Join[{pf[q][[1]]}, magSym syms]]],
          spatKeys];
    (* spatial vectors LAST, off the finished components: whatever treatment the parent got — a
       polyFrameSpec external, a full unit loop, or a spatial unit loop — the spatial vector is that
       same component list with the temporal slot zeroed, sharing the parent's ntU$ symbols and its
       unit group. No new symbol, no new group. *)
    With[{done = Join[KeyDrop[pf, spatKeys], nf, nfS]},
      {Join[done, Association[(# -> ReplacePart[done[First[#]], 1 -> 0]) & /@ svKeys]],
       defs, groups}]];

(* Whether `frame` matches the unit-loop spec's assumption: every momentum is either an external
   depending only on `pSym`, or a loop whose every component is `dir·magSym` (proportional to the
   single magnitude). A finite-T frame breaks this — the external carries an independent temporal
   component p0 and the loop an independent l0 (neither ∝ magSym) — so we must fall back to the
   general polyFrameSpec. *)

unitLoopOkQ[frame_, pSym_, magSym_] := !ntEnvFlag["NT_NO_UNIT_GROUPS"] && AllTrue[
    Values[frame],
    Function[comps,
      Module[{cc = PowerExpand[comps]},
        SubsetQ[{pSym}, Variables[cc]] || fullLoopQ[cc, magSym]
      ]]];

numericComponents[env_, frame_, symDefs_, unitGroups_ : {}] := Module[
    {compExpr, usyms, nsym, mpcpp, compCpp, varFill, vfill, idx, units},
    (* 4 components per momentum Base (polynomial in the user symbols). *)
    compExpr = Association @ KeyValueMap[#2["Base"] -> resolveComponents[#1, frame]&, env];
(* reject only fractional powers of NON-numeric bases (e.g. Sqrt[1-cos1^2] not rewritten to a sin
   symbol); a numeric irrational coefficient like Sqrt[3]/2 from the 120-degree external frame is a
   perfectly good polynomial coefficient. *)
    With[{bad = DeleteDuplicates @ Cases[Values[compExpr], Power[b_, _Rational] /; !NumericQ[b], Infinity]},
      If[bad =!= {},
        Print["numericComponents: non-polynomial components (fractional power of a symbol remains): ", bad];
        Abort[]]];
    usyms = Sort @ DeleteDuplicates @ Flatten[Variables /@ Values[compExpr]];
    nsym = Length[usyms];
(* Coefficients may be complex (a silver-blaze frame component pi T - I muq has coefficient -I on
   muq): emit Re and Im separately, as dscV does -- cppNum of a Complex is Wolfram's Complex(a,b),
   which is not C++. Fixed 2026-09-24. *)
(* Emit through the generator's `env` (a LorentzEnv bound to nsym) — the sole construction path for
   MPoly now that the bare-nsym factories are private. *)
    mpcpp[e_] := Module[{rules = CoefficientRules[e, usyms]},
        If[rules === {},
          "env.zero()",
          "(" <> StringRiffle[("env.mono({" <> StringRiffle[ToString /@ #[[1]], ","] <> "},Cx{" <> cppNum[Re[#[[2]]]] <> "," <> cppNum[Im[#[[2]]]] <> "})")& /@ rules, " + "] <> ")"
        ]];
    compCpp = Association @ KeyValueMap[#1 -> (mpcpp /@ #2)&, compExpr];
    vfill[s_] := If[KeyExistsQ[symDefs, s],
        cppFlat[symDefs[s]],
        SymbolName[s]];
    varFill = vfill /@ usyms;(* indexed by MPoly var id (0-based) *)
(* Unit-constraint groups (ΣUμ²=1) as MPoly var-index lists, so the C++ reduce_units collapses the
   bare-loop denominator to the monomial l1² and the U·U projector factors. The caller supplies the
   groups as symbol lists (the loop's unit-direction components); drop any symbol not in usyms (a
   component that vanished) and any group with <2 surviving symbols. *)
    idx = AssociationThread[usyms -> Range[Length[usyms]] - 1];
    units =
      DeleteCases[
        Function[g,
            Lookup[idx, Select[g, KeyExistsQ[idx, #]&]]
          ] /@ unitGroups,
        _ ? (Length[#] < 2&)];
    <|
      "nsym" -> nsym,
      "usyms" -> usyms,
      "compCpp" -> compCpp,
      "varFill" -> varFill,
(* default -1 (-> comp size 0) for a purely scalar integrand whose 4-vector component env is
   empty: Max[{}] is -Infinity, which would leak into the C++ as the comp() vector size. The loop
   momentum of such a flow (e.g. a bosonic meson-potential tadpole) is carried by usyms/units. *)
      "symNamesCpp" -> (("(" <> # <> ")")& /@ varFill),
      "maxBase" -> Max[Append[#["Base"]& /@ Values[env], -1]],
      "units" -> units
    |>];
