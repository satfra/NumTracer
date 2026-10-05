(* CodegenFrames.m — momentum frames to polynomial component tables. Rewrites a kinematic frame into
   components polynomial in fresh symbols (unit-loop, mixed, or general polyFrameSpec), and builds the
   numeric backend's component table (numericComponents). Loaded by NumTracer.m via ntLoadPart, in
   NumTracer`Private`. *)

(* SPATIAL VECTORS (ntSpatialVec[q], from FormTracer's vecs) are not independent momenta: their
   components are the parent's with the temporal slot zeroed. Minting them their own unit group would
   be numerically correct but hides the identity from reduce_units, so cross terms never collapse
   against the parent's Σ U² = 1. They are held out of the classification and derived from the
   finished parent components. Select, not Complement: the key ORDER fixes the ntU$ numbering. *)
spatialVecKeysOf[frame_] :=
  Select[Keys[frame], MatchQ[#, _ntSpatialVec] && KeyExistsQ[frame, First[#]] &];

(* Unit-loop spec, the compact parametrisation for vacuum frames. Writing the loop as
   MAGNITUDE × UNIT-DIRECTION, comp_μ = l1 · Uμ with degree-1 symbols Uμ (fill computes Uμ = dirμ),
   makes l·p_i degree-1 in Uμ, and the bare-loop denominator l1²·ΣUμ² collapses to the monomial l1²
   under the unit constraint ΣUμ² = 1 (returned as a unit `group`), so its 1/l² atom cancels. The
   general polyFrameSpec instead expands the loop into angle products whose degrees compound through
   a projector chain. Externals depend only on `pSym`; the loop is the momentum carrying `magSym`. *)
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

(* General frame -> {polyFrame, defs}: replace the trig and radical sub-expressions in the frame
   components with fresh symbols so the components become polynomials. Introduces a symbol for each
   loop polar factor Sqrt[1-cos^2] (a sin), each bare-angle Cos/Sin, and each remaining radical;
   `defs` maps each new symbol back to its closed form so the kernel can fill it. An already
   polynomial "Components" association passes through unchanged. The fallback when neither unit-loop
   spec applies. *)
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
    (* GENERAL radicals of composite arguments (e.g. Sqrt[S0^2 (1 - S1 Sin[SPhi])] in the (S0,S1,SPhi)
       frame), which numericComponents would otherwise reject. One fresh symbol per distinct
       (base, exponent) pair, so negative half-powers become plain symbols too. ReplaceAll, not
       ReplaceRepeated: it rewrites the OUTERMOST radical only, so a nested radical stays inside the
       recorded closed form and is emitted as part of that one C++ expression. *)
    Module[{n = 1, radSeen = <||>},
      polyFrame =
        polyFrame /. (pw : Power[b_, _Rational] /; ! NumericQ[b]) :>
          (If[! KeyExistsQ[radSeen, pw],
             With[{new = Symbol["ntRad$" <> ToString[n++]]},
               radSeen[pw] = new;
               (* record the closed form in the raw runtime coordinates: ntSin$/ntCos$/ntSinA$ are
                  fill slots, not C++ variables, and cppFlat emits a def verbatim *)
               defs[new] = pw //. Normal[defs]]];
           radSeen[pw])];
    polyFrame = Association @ Thread[Keys[frame] -> polyFrame];
    {polyFrame, defs}];

(* a frame component that is magSym times a magSym-free coefficient *)
magPropQ[c_, magSym_] := Simplify[c - Coefficient[c, magSym] magSym] === 0;
(* all four components of an (already PowerExpand-ed) frame vector are magSym-proportional *)
fullLoopQ[cc_, magSym_] := AllTrue[Range[4], magPropQ[cc[[#]], magSym]&];

(* SPATIAL unit loop (finite T): the temporal component (slot 1 ↔ C++ component 0) is an independent
   magSym-free coordinate (l0 / a Matsubara frequency), while the SPATIAL components are all ∝ magSym
   with at least one nonzero direction. O(3) is intact at T > 0, so the spatial direction is a unit
   3-vector and Σ U² = 1 holds over the spatial components; the bare-loop denominator reduces to the
   two-term l0² + l1². Assumes the frame parametrises the spatial part polar-wise (unit-norm dirs). *)
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

(* MIXED unit-loop spec: loop momenta as magnitude × unit-direction, externals (which may carry shape
   coordinates, e.g. the (S0,S1,SPhi) frames) through polyFrameSpec's radical-minting. The loop is
   frame-independent (∝ l1 in every frame), so it keeps the unit-direction treatment regardless of
   the externals; sending the whole frame through polyFrameSpec makes the traces far larger. Keeping
   polyFrameSpec's cos/sin symbols with {cos, ntSin$cos} pair groups instead is WORSE still (sin² →
   1 − cos² expands monomials binomially); the gain comes from the loop entering at degree 1.
   Loop tags with the same direction coefficients (e.g. l1 and its fermionic partner lf1 = l1 +
   (pi T, 0, 0, 0)) share one ntU$ set and one unit group, so their monomials merge and reduce
   against the same Σ U² = 1; `dirSyms` memoises by (kind, directions).
   NT_NO_UNIT_GROUPS sends the whole frame through polyFrameSpec (grading control and rollback). *)
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

(* ---- numeric (matrix-product) backend: component table + user symbols ----------------------------
   The numeric backend evaluates scalar products from each momentum's 4 COMPONENTS, so it needs only:
   the polynomial variables (the free user symbols), each fundamental momentum's 4 components as
   Poly-builder C++, and the C++ fill formula for each symbol (a kernel argument, or a derived
   symbol like sin1 = sqrt(1-cos1^2)). Composite momenta resolve by component arithmetic via
   resolveComponents. *)
numericComponents::nonpoly = "Non-polynomial momentum components (a fractional power of a symbol remains): `1`";

numericComponents[env_, frame_, symDefs_, unitGroups_ : {}] := Module[
    {compExpr, usyms, nsym, mpcpp, compCpp, varFill, vfill, idx, units},
    (* 4 components per momentum Base (polynomial in the user symbols). *)
    compExpr = Association @ KeyValueMap[#2["Base"] -> resolveComponents[#1, frame]&, env];
    (* reject only fractional powers of NON-numeric bases (e.g. Sqrt[1-cos1^2] not rewritten to a
       sin symbol); a numeric irrational coefficient like Sqrt[3]/2 is a valid coefficient *)
    With[{bad = DeleteDuplicates @ Cases[Values[compExpr], Power[b_, _Rational] /; !NumericQ[b], Infinity]},
      If[bad =!= {},
        Message[numericComponents::nonpoly, bad];
        Abort[]]];
    usyms = Sort @ DeleteDuplicates @ Flatten[Variables /@ Values[compExpr]];
    nsym = Length[usyms];
    (* Coefficients may be complex (a silver-blaze component pi T - I muq), so emit Re and Im
       separately: cppNum of a Complex is Wolfram's Complex(a,b), which is not C++. A Poly is built
       only through the generator's `frame` (a Frame bound to the symbol list). *)
    mpcpp[e_] := Module[{rules = CoefficientRules[e, usyms]},
        If[rules === {},
          "frame.zero()",
          "(" <> StringRiffle[("frame.mono({" <> StringRiffle[ToString /@ #[[1]], ","] <> "},Cx{" <> cppNum[Re[#[[2]]]] <> "," <> cppNum[Im[#[[2]]]] <> "})")& /@ rules, " + "] <> ")"
        ]];
    compCpp = Association @ KeyValueMap[#1 -> (mpcpp /@ #2)&, compExpr];
    vfill[s_] := If[KeyExistsQ[symDefs, s],
        cppFlat[symDefs[s]],
        SymbolName[s]];
    (* indexed by Poly var id (0-based) *)
    varFill = vfill /@ usyms;
    (* Unit-constraint groups (ΣUμ²=1) as Poly var-index lists, so the C++ reduce_units collapses
       the bare-loop denominator to the monomial l1². Drop symbols not in usyms (a vanished
       component) and groups with < 2 surviving symbols. *)
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
      "symNamesCpp" -> (("(" <> # <> ")")& /@ varFill),
      (* default -1 (comp size 0) for a purely scalar integrand with an empty component env:
         Max[{}] is -Infinity, which would leak into the C++ as the comp() vector size *)
      "maxBase" -> Max[Append[#["Base"]& /@ Values[env], -1]],
      "units" -> units,
      (* a finite-density frame (p0 - I muq): projector denominators may be complex *)
      "complex" -> !FreeQ[Values[compExpr], Complex]
    |>];
