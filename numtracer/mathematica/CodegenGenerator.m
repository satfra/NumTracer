(* ---- emit the NUMERIC generator program. Builds each net's DiracNet chain + Lorentz/projector
        NetVal in C++, contracts them numerically (4×4 matrix products), folds colour per group, and
        PRINTS the committed straight-line kernel header. No reduce/rebase/ibp/sp-kinematics. *)

(* ---- STAGE: the global sub-term dedup JOIN ----------------------------------------------------
   Input: the six already-interned id columns (one ragged list per net) and the five pools they index.
   Output: the distinct traces, the per-net fold entries that reference them, and the counts the
   caching decision needs. See the "GLOBAL SUB-TERM DEDUP" note above for WHY this exists.

   THE FOUR INVARIANTS THAT MAKE THE OUTPUT BYTE-IDENTICAL. Every one of them fixes an ORDER, and the
   emitted kernel's slot numbering follows that order, so breaking one produces a kernel that is
   still correct and still passes every oracle while differing byte-for-byte from the committed one:

     I1  The pools arrive in ntMkIntern's FIRST-APPEARANCE order, and this stage never reorders them.
     I2  The mixed-radix pack reads its digits least-significant-last —
           traceKeyPacked = ((ds*nLs + ls)*nDc + dc)*nDl + dl,  and then *nDr + dr for the channel —
         and the decode below reverses exactly that order. Radices are Max[1, Length[pool]]: an empty
         column would otherwise give a zero radix and a degenerate pack.
     I3  PositionIndex groups in first-appearance order and ReverseSort on an Association is STABLE.
         Do not "simplify" either to GatherBy / SortBy — both would regroup.
     I4  Under NT_GEN_NO_DEDUP the grouping id becomes a flat global occurrence index instead of the
         packed key, so nothing merges. The substitution and the packOfDistinct read-back are ONE
         unit: two committed reference kernels are generated this way and graded against the
         deduped ones (tests/refshim/compare_lambda3d_small.cpp, compare_zaaqbq1_small.cpp). *)

(* One net's sub-terms -> {trace keys, dress ids, summed scalars}, one entry per distinct
   (trace, dress channel) key in first-appearance order (I3: PositionIndex), zero sums dropped.
   A packed numeric scalar column takes the vectorised path: singletons are copied, pairs are added
   elementwise (identical to Total of two), larger groups use Total. Anything else (a symbolic
   scalar) keeps the per-group path, which keeps every sum not provably zero. *)
ntMergeNetTerms[ck_, tks_, drs_, scs_] := Module[{g = Values[PositionIndex[ck]], f, sums, len, keep},
  f = g[[All, 1]];
  If[Developer`PackedArrayQ[scs],
    len = Length /@ g;
    sums = scs[[f]];
    With[{p2 = Flatten @ Position[len, 2, {1}]},
      If[p2 =!= {}, sums[[p2]] = scs[[g[[p2, 1]]]] + scs[[g[[p2, 2]]]]]];
    With[{pk = Flatten @ Position[UnitStep[len - 3], 1, {1}]},
      If[pk =!= {}, sums[[pk]] = Total[scs[[#]]] & /@ g[[pk]]]];
    (* numeric, not structural: a packed machine sum cancels to 0. + 0. I, which =!= 0 would keep *)
    keep = Unitize[Abs[sums]];
    {Pick[tks[[f]], keep, 1], Pick[drs[[f]], keep, 1], Pick[sums, keep, 1]},
    sums = Total[scs[[#]]] & /@ g;
    keep = If[TrueQ[# == 0], 0, 1] & /@ sums;
    {Pick[tks[[f]], keep, 1], Pick[drs[[f]], keep, 1], Pick[sums, keep, 1]}]];

ntGenDedupJoin::bigkey = "The packed sub-term key range `1` exceeds 2^62; the dedup join stays exact but runs on unpacked bignum columns (slow). Consider re-ranking the key columns to dense ids.";

ntGenDedupJoin[diracNetIds_, lorNetIds_, subScalars_, dressChainIds_, slotTupleIds_, dressMonoIds_,
               diracNetPool_, lorNetPool_, dressChainPool_, slotTuplePool_, dressMonoPool_,
               hasDressed_, noDedup_] :=
  Module[{lens, dsI, uDs, lsI, uLs, dcI, uDc, dlI, uDl, drI, uDr,
          nLs, nDc, nDl, nDr, traceKeyPacked, traceDressKeyPacked, traceKeyFlat, distinctTraceKeys,
          subKeysLen, merged, foldKeys, foldLens, netTraceRows, netDressRows, netScalarRows,
          refCount, distinctSubs, subIdxOf, nSub, nReused},
(* The five key columns arrived already interned (see the expansion above), so the join opens with
   the ids and the pools in hand: no DeleteDuplicates, no Lookup, nothing hashed here at all. This is
   what the stage used to spend most of its time on — 5 columns x 12.7 M elements x 2 passes. *)
    lens = Length /@ diracNetIds;

    {dsI, uDs} = {diracNetIds, diracNetPool};
    {lsI, uLs} = {lorNetIds, lorNetPool};
    (* a non-dressed flow has no chain/slot column: the scalar 0 broadcasts in the listable pack *)
    {dcI, uDc} = If[hasDressed, {dressChainIds, dressChainPool}, {0, {""}}];
    {dlI, uDl} = If[hasDressed, {slotTupleIds, slotTuplePool},   {0, {""}}];
    {drI, uDr} = {dressMonoIds, dressMonoPool};
    {nLs, nDc, nDl, nDr} = Max[1, Length[#]]& /@ {uLs, uDc, uDl, uDr};   (* I2 *)
    (* past 2^63 the packed keys become bignums: still exact, but the columns unpack and slow down *)
    If[Max[1, Length[uDs]] nLs nDc nDl nDr >= 2^62,
      Message[ntGenDedupJoin::bigkey, Max[1, Length[uDs]] nLs nDc nDl nDr]];

(* I2: mixed-radix pack of the traceKey, then of the (traceKey, dressChannel) pair the per-net merge
   groups on. Listable arithmetic over the ragged integer columns — no Map. *)
    traceKeyPacked      = ((dsI * nLs + lsI) * nDc + dcI) * nDl + dlI;
    traceDressKeyPacked = traceKeyPacked * nDr + drI;
    subKeysLen = lens;
    If[noDedup,   (* I4 *)
      Module[{tot = Total[lens]},
        traceKeyFlat = Flatten[traceKeyPacked, 1];
        traceKeyPacked = TakeList[Range[0, tot - 1], lens];
        traceDressKeyPacked = traceKeyPacked],
      traceKeyFlat = None];

(* Per net: merge the sub-terms that share a trace AND a dress channel, summing scalars; drop the zero
   sums. Do this BEFORE counting references — a (trace, channel) occurring twice inside ONE net collapses
   to a single reference here, so its raw occurrence count would overstate its reuse, and a channel whose
   scalars cancel is not referenced at all and must not be contracted. Each net yields three columns
   {trace key, dress id, summed scalar} (ntMergeNetTerms). *)
    merged = MapThread[ntMergeNetTerms, {traceDressKeyPacked, traceKeyPacked, drI, subScalars}];
    foldKeys = Join @@ merged[[All, 1]];
    foldLens = Length /@ merged[[All, 1]];
(* I3: order the distinct traces by DESCENDING reference count, so a memory-capped run still caches
   the traces that repay caching most, and the singletons (refCount 1 — computing one costs the same
   whether or not it is cached, so caching it is pure RAM for no saving) land at the end where the
   default cap excludes them. ReverseSort on an Association is stable, which is what keeps ties in
   first-appearance order. *)
    refCount = Counts[foldKeys];
    distinctSubs = Keys[ReverseSort[refCount]];
    subIdxOf = AssociationThread[distinctSubs -> Range[0, Length[distinctSubs] - 1]];
    nSub = Length[distinctSubs];
    nReused = Total[UnitStep[Values[refCount] - 2]];(* #refCount >= 2 == the length of the sorted prefix worth caching *)
(* I2/I4: decode the surviving packed keys back into their four pool entries, listable over all nSub
   keys at once. Under NT_GEN_NO_DEDUP the grouping id is an occurrence index, so the pack is read off
   the flat column. *)
    distinctTraceKeys = If[traceKeyFlat === None, distinctSubs, traceKeyFlat[[distinctSubs + 1]]];
    distinctSubs =
      Module[{a = distinctTraceKeys, diracIdx, lorIdx, chainIdx, slotIdx},
        slotIdx  = Mod[a, nDl]; a = Quotient[a, nDl];
        chainIdx = Mod[a, nDc]; a = Quotient[a, nDc];
        lorIdx   = Mod[a, nLs]; diracIdx = Quotient[a, nLs];
        Transpose[{uDs[[diracIdx + 1]], uLs[[lorIdx + 1]], uDc[[chainIdx + 1]], uDl[[slotIdx + 1]]}]];
(* the fold entries, per net: distinct-trace index, the dress-atom multiset itself (downstream builds
   the DMono table from it) and the scalar *)
    netTraceRows  = TakeList[Lookup[subIdxOf, foldKeys], foldLens];
    netDressRows  = TakeList[uDr[[(Join @@ merged[[All, 2]]) + 1]], foldLens];
    netScalarRows = merged[[All, 3]];

    ntStageResult["ntGenDedupJoin",
      {"subKeysLen", "netTraceRows", "netDressRows", "netScalarRows", "distinctSubs", "nSub", "nReused"},
      <|"subKeysLen" -> subKeysLen, "netTraceRows" -> netTraceRows, "netDressRows" -> netDressRows,
        "netScalarRows" -> netScalarRows,
        "distinctSubs" -> distinctSubs, "nSub" -> nSub,
        "nReused" -> nReused|>]];

(* ---- emitNumericGenerator: STAGE MAP ----------------------------------------------------------
   Emits the build-time generator PROGRAM — a main TU, N net-builder unit TUs, and the decl header
   they share. Returns {pre, units, decl, main}. The program it emits is what actually contracts the
   traces and PRINTS the committed straight-line kernel header; nothing here contracts anything.

   The stages, in order. Each consumes the one before it; the numbers are the shape of the data, not
   an aspiration:

     1. SUB-TERM EXPANSION + INTERNING. Walk every net's cores, expand the dressed slot options into
        their Cartesian product, and intern the five key columns (Dirac net, Lorentz net, dressed
        chain, slot-option tuple, dressing-atom multiset) through ntMkIntern. Output: six ragged
        integer columns, one entry per (net, sub-term), plus the five pools they index.
     2. NET CSE. Net-builder strings that recur across sub-terms become shared `lc<k>()` / `dc<k>()`
        accessor functions, and the pools are rewritten to call them.
     3. DEDUP JOIN (ntGenDedupJoin, above). Merge sub-terms sharing a trace AND a dress channel,
        count references, and order the distinct traces by descending reference count. Output: the
        distinct-trace table and the per-net fold entries that index it. This is the stage the whole
        design exists for — a dense flow has 5-7x more sub-terms than distinct traces.
     4. UNIT TABLES. Chunk the distinct-trace tables (ntChunkDefs) and LPT bin-pack every definition
        into the -O0 unit TUs, so the net-builder code compiles in parallel.
     5. main() DATA TABLES. Emit the per-net index/scalar/dressing tables with two-level hash-consing
        (distinct ROWS, then distinct VALUES), each as a top-level chunk function collected into
        `tableBag` and prepended afterwards.
     6. PHASE A / PHASE B / EMISSION. Emit the calls that contract each distinct trace (phase A),
        fold each net and drain each group (phase B), lower to SSA, and print the kernel header.

   WHAT MAKES THE OUTPUT REPRODUCIBLE. Every stage boundary above fixes an ORDER, and the emitted
   kernel's slot numbering follows it, so a reordering here yields a kernel that is still correct and
   still passes every oracle while differing byte-for-byte from the committed one. The four specific
   invariants are stated at ntGenDedupJoin (I1-I4); the fifth lives at stage 4: `allDefs` is Joined in
   a fixed order and `Ordering` breaks ties by index, so that Join order decides which definition
   lands in which unit TU.

   The `Tuples` in stage 1 runs LAST-SLOT-FASTEST, and stages 1, 3 and 5 all rely on that alignment.
   It is the one convention worth having in mind while reading any of them. *)
(* `mIdx` is the MPoly var index (0-based) of the Matsubara frequency, or -1 for "not a finite-T
   flow / unknown". When it is >= 0 the generator proves Matsubara evenness while it contracts (see
   the ntMEven thread below) and emits the verdict as a constant in the traces header. *)
emitNumericGenerator[invNets_, invRest_, colourNets_, groups_, ncomp_, nsInner_, fillArgSig_, kns_:"numtracer_kernels", complexQ_:False, realOnlyG_ : {}, crossCSE_:False, mIdx_:-1] :=
  Module[{nNet = Length[invNets], nGrp = Length[groups], nsym = ncomp["nsym"], maxBase = ncomp["maxBase"], varFill = ncomp["varFill"], symNames = ncomp["symNamesCpp"], compCpp = ncomp["compCpp"], unitG = ncomp["units"], str, tmpl, pre, unitPre, unitInc, nUnits, units, decl, diracNetStrs, lorentzNetStrs, subScalars, dressChains, dressSlotOpts, hasDressed, allDefs, main, compInit, cseDefs, cseDecls, chunkDecls, ntNoDedup, subKeysLen, netTraceRows, netDressRows, netScalarRows, dsInt, dsGet, lsInt, lsGet, dcInt, dcGet, dlInt, dlGet, drInt, drGet, scCache, slotCombo, dPoolCse, lPoolCse, distinctSubs, nSub, nReused, sdnDefs, slnDefs, sdnCDecl, slnCDecl, tableBag, emitBigTable,
    chpDefs, chpCDecl, optpDefs, optpCDecl, sdchrDefs, sdchrCDecl, sdslrDefs, sdslrCDecl, dressAtomIds, colMainDecls, colChunkDefs},
    str[x_] := ToString[x];
(* The four big index/scalar tables and the colour/group tables are emitted as TOP-LEVEL chunk
   functions (ntBigTableFns) that main() then calls, because inline they make main() large enough to
   break the compiler. Their DEFINITIONS therefore have to be prepended to main() after the body is
   built — so the body's construction collects them as a side effect, into this bag.

   The bag (not a string being rebuilt with <>) makes the accumulation explicit and O(1) per append.
   The ORDER is still fixed by Mathematica's left-to-right evaluation of the StringJoin arguments
   below, and that order is load-bearing: it decides the order the table functions appear in the
   generator source. Byte-identical to the string accumulation it replaces, for exactly that reason. *)
    tableBag = Internal`Bag[];
(* Emit one big table: stuff its definition into the bag, and return the `<lhs> = <call>;` line that
   goes in main(). Collapses seven copies of the same three-line With/mutate/return shape. *)
    emitBigTable[nm_String, ret_String, rows_List, lhs_String] :=
      With[{t = ntBigTableFns[nm, ret, rows]},
        Internal`StuffBag[tableBag, t[[1]]];
        "  " <> lhs <> " = " <> t[[2]] <> ";\n"];
(* DRESSED nets (symbolic dressing collection): a core may be ntDressedCore[chainStr, slotsStr]
   (a numerator structure-sum kept eager). LEVER (b): the trace table is PLAIN MPoly for both paths — a
   dressed sub-term's structural trace (dressing stripped) contracts via numeric_value_dressed_netval_mp,
   and its dressing rides the per-sub-term scalar (dsc, numeric) + monomial (sdr, atom ids). The phase-B
   fold (fold_groups_streaming_dressed) assembles those into the per-net DPoly, and the dressings ride the
   env as kind-2 `dress` leaves filled by fm.dress. So combinations that share a concrete structure but
   differ only in dressing collapse to ONE trace; the non-dressed path is byte-identical. *)
    hasDressed = !FreeQ[invNets, _ntDressedCore];
    (* shared wrapper templates so the net-builder strings compile. *)
    tmpl = "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal tproj(){ return projT(Mu,Nu,Lb,Inv); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal lproj(){ return projL(Mu,Nu,Lb,Inv); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int InvS> NetVal mproj(){ return projM(Mu,Nu,Lb,InvS); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int Inv,int InvS> NetVal eproj(){ return projE(Mu,Nu,Lb,Inv,InvS); }\n" <> "template<int Mu,int Nu> NetVal lmetric(){ return met(Mu,Nu); }\n" <> "template<int Lbl,int Base,int Mask> NetVal lvec(){ return vec(Lbl,Base); }\n" <> "template<int A,int B,int C,int D> NetVal leps(){ return epsilon(A,B,C,D); }\n" <> "inline NetVal konst(double c){ return NetVal{PTerm{Cx{c,0}, {}}}; }\n" <> "template<class L> struct litco;\n" <> "template<numtracer::Cx C> struct litco<numtracer::Lit<C>>{ static constexpr numtracer::Cx v=C; };\n" <> "template<class L> NetVal sc(NetVal x){ return scale(litco<L>::v, std::move(x)); }\n";
(* per net: a colour group is a SUM of sub-terms. invNets[i] = {core_b…} (each a DiracNet literal
   for a gamma branch, or a Lorentz NetVal for a gamma-free branch); invRest[i] = {{rest_b,scal_b}…}
   parallel. Build per net the parallel lists of {DiracNet builder, NetVal builder} (a gamma-free
   branch → empty DiracNet + the whole net as the rest) plus the sub-term scalars; the generator
   sums mp[i] = Σ_b scal_b · numeric_value_netval(dn[i][b], ln[i][b]). *)
(* Each branch yields a LIST of {ds, ls, scal, dc, dl, dr} sub-terms (usually length 1). A DRESSED branch
   expands the Cartesian product of its chain's slot options (Tuples) into ONE structural sub-term per
   combination. LEVER (b): each slot option is now a TRIPLE {structStr, num, dr} (dressing-free DSlotOpt +
   numeric Cx + dress-atom ids). For each combination we keep the STRUCTURE (dl = the list of structStr,
   one per slot) for the trace key, fold the numeric part (∏ num) INTO the sub-term scalar, and carry the
   dressing-atom multiset (dr = ⋃ dress) as the DPoly key. So combinations that share a concrete structure
   but differ only in dressing collapse to ONE plain-MPoly trace (the 6.2× dedup), and the trace table
   loses its dressing dimension entirely. A slot's option list is nv[[2]][[k]]; Tuples over them gives every
   combination. The non-dressed branches carry an empty dress key, so non-dressed flows stay byte-identical. *)
(* Timed: this is where the dressed slot options are expanded into their Cartesian product, so it is
   the first place in the emitter that touches every SUB-TERM individually (millions of them on a
   dressed flow) rather than every net. Separate from the dedup join below because the work is
   different in kind — construction, not hashing. *)
(* INTEGER KEY COLUMNS. The five key columns (dirac net, lorentz net, chain, structural option tuple,
   dress multiset) used to be carried as STRINGS and lists-of-strings all the way to the dedup join,
   which then hashed each of them twice more (DeleteDuplicates + Lookup, per column). On ZAAqbq2 that
   is 12.7 M elements x 5 columns x 3 passes. Intern them HERE instead, once per distinct value, and
   carry ids: the join's `internCol` disappears and its arithmetic is over packed integers.
   Byte-identity: `ntMkIntern` assigns ids in first-appearance order, and the walk order here
   (nets -> cores -> combinations, `Tuples` last-slot-fastest) is exactly the order the flattened
   column had, so the id -> value map is the same permutation `DeleteDuplicates` produced. Each
   interner is called ONCE per core-instance (or once per distinct slot-option set, via slotCombo)
   rather than once per sub-term, which is the 33x. *)
    {dsInt, dsGet} = ntMkIntern[];
    {lsInt, lsGet} = ntMkIntern[];
    {dcInt, dcGet} = ntMkIntern[];
    {dlInt, dlGet} = ntMkIntern[];
    {drInt, drGet} = ntMkIntern[];
(* The Cartesian expansion of one core's slot options depends on `slotOpts` ALONE — not on the chain,
   the Lorentz rest or the branch scalar. So compute it once per distinct option set and reuse: on
   ZAAqbq2 that is ~36 k evaluations instead of ~1.27 M. Returns {dlIds, drIds, nums, n}, all four
   from one call so the columns can never desync — `Tuples` and `Flatten[Outer[...]]` both vary the
   last slot fastest, and that alignment between a structure tuple and its numeric coefficient is
   load-bearing. *)
    scCache = <||>;
    slotCombo[slotOpts_] :=
      Lookup[scCache, Key[slotOpts],
        scCache[slotOpts] =
          Module[{structs, nums, dress, cs, n},
            structs = #[[All, 1]]& /@ slotOpts;
            nums    = #[[All, 2]]& /@ slotOpts;
            dress   = #[[All, 3]]& /@ slotOpts;
            cs = Tuples[structs];
            n  = Length[cs];
            {Developer`ToPackedArray[dlInt /@ cs],
(* by far the common case: no slot option on this core carries a dressing atom, so every
   combination's union is empty. Checked once per distinct option set. *)
             If[AllTrue[dress, AllTrue[#, # === {}&]&],
               ConstantArray[drInt[{}], n],
               Developer`ToPackedArray[drInt /@ (Sort[Catenate[#]]& /@ Tuples[dress])]],
             Flatten[Outer[Times, Sequence @@ nums]],
             n}]];
    With[{ntT = First @ AbsoluteTiming[
    {diracNetStrs, lorentzNetStrs, subScalars, dressChains, dressSlotOpts, dressAtomIds} =
      Transpose @
        MapThread[
          Function[{cores, rss},
            If[cores === {},
              {{}, {}, {}, {}, {}, {}},
(* COLUMN-ORIENTED, not row-oriented. The obvious spelling builds one six-element ROW per sub-term
   and Transposes at the end — but a dressed flow has millions of sub-terms (12.7 M on ZAAqbq2), and
   materialising a row for each, with three separate `#[[k]]& /@ combo` passes inside it, measured
   63 s. Every column is either a constant repeated n times, a single `Tuples`, or a single `Outer`
   product, so build the six columns directly and concatenate them per net; nothing is evaluated per
   sub-term at Mathematica level. `Join @@@ Transpose[...]` regroups the per-core column tuples into
   six whole-net columns.
   ORDER is what makes this exact: `Tuples` varies the LAST slot fastest and `Flatten[Outer[...]]`
   does the same, so the numeric product lines up element-for-element with the structure tuple it
   belongs to. Verified byte-identical on the dressed fixtures. *)
              Join @@@ Transpose @
                  MapThread[
                    Function[{nv, rv, scal},
                      Module[{lsStr = If[rv === "", "NetVal{}", rv]},
                        Which[
                          MatchQ[nv, _ntDressedCore],(* dressed numerator: expand slot options → structural sub-terms *)
                            With[{chain = nv[[1]], slotOpts = nv[[2]]},
                              If[slotOpts === {},
                                {{dsInt["DiracNet{}"]}, {lsInt[lsStr]}, ntPackCx[{scal}], {dcInt[chain]}, {dlInt[{}]}, {drInt[{}]}},
(* dl = this combination's STRUCTURAL option-string LIST (one dressing-free structStr per chain
   slot), interned so the table emitter still pools the distinct structures; the numeric Cx folds
   into the sub-term scalar and the dress ids become the DPoly key. *)
                                Module[{tb = slotCombo[slotOpts], n},
                                  n = tb[[4]];
                                  {ConstantArray[dsInt["DiracNet{}"], n],
                                   ConstantArray[lsInt[lsStr], n],
                                   ntPackCx[scal * tb[[3]]],
                                   ConstantArray[dcInt[chain], n],
                                   tb[[1]],
                                   tb[[2]]}]]],
                          StringStartsQ[nv, "DiracNet"],(* gamma branch: DiracNet + projector rest *)
                            {{dsInt[nv]}, {lsInt[lsStr]}, ntPackCx[{scal}], {dcInt["std::vector<DChainTok>{}"]}, {dlInt[{}]}, {drInt[{}]}},
                          True,(* gamma-free branch: whole net is the rest *)
                            {{dsInt["DiracNet{}"]}, {lsInt[nv]}, ntPackCx[{scal}], {dcInt["std::vector<DChainTok>{}"]}, {dlInt[{}]}, {drInt[{}]}}]]],
                    {cores, rss[[All, 1]], rss[[All, 2]]}]]],
          {invNets, invRest}];]},
      ntLog["[prof] sub-term expansion: ", ntT, " s"]];
(* ---- colour-net table: chunk DEFINITIONS on the parallel -O0 units, assembler in the main TU ----
   The distinct colour nets are one `SUNNet{sun3.T(..), ..}` constructor-call literal each, and on a
   flow with a large colour graph the table dwarfs everything else in the main TU (measured: 6.28 MB
   of a 9.94 MB TU, 15226 distinct nets) — and the main TU is the ONE -O1 compile that cannot be
   parallelised, so the table sat squarely on the compile critical path. The chunk functions have
   exactly the net-builder shape, so they now ride the same LPT-packed -O0 units as every other
   builder (colChunkDefs -> allDefs below); only forward decls + the tiny assembler stay in the main
   TU. External (non-static) linkage is what makes the cross-TU split work. Chunking rationale and
   history (clang -O1 >1h / -O0 stack overflow on the braced-init inside main) in git. Element order
   is preserved, so the assembled vector is element-for-element what the braced-init produced. *)
    {colMainDecls, colChunkDefs} =
      With[{uCol = DeleteDuplicates[colourNets]},
        If[uCol === {},
          {"static std::vector<SUNNet> ntColNets(){ return {}; }\n", {}},
          Module[{envDecl, chunks},
            envDecl = StringJoin["SUNEnv sun" <> # <> "(" <> # <> "); "& /@ DeleteDuplicates @ Flatten @ StringCases[uCol, "sun" ~~ r : DigitCharacter.. ~~ "." :> r]];
            chunks = ntSplitByChars[uCol, 2];
            {StringJoin[
               MapIndexed["void ntColNets_c" <> ToString[#2[[1]] - 1] <> "(std::vector<SUNNet>& o);\n"&, chunks],
               "static std::vector<SUNNet> ntColNets(){ std::vector<SUNNet> o; o.reserve(" <> ToString[Length[uCol]] <> "); " <> StringJoin[Table["ntColNets_c" <> ToString[k - 1] <> "(o); ", {k, Length[chunks]}]] <> "return o; }\n"],
             MapIndexed["void ntColNets_c" <> ToString[#2[[1]] - 1] <> "(std::vector<SUNNet>& o){ " <> envDecl <> StringJoin[("o.push_back(" <> # <> "); ")& /@ #1] <> "}\n"&, chunks]}]]];
    pre =
      StringJoin[
        "// GENERATED by MakeNTKernel — do not edit. Numeric matrix-product tensor traces.\n",
(* one umbrella header pulls the whole engine API (network/dirac/gen/sun_net/numeric_*/core) for
   this single generator TU. The parallel -O0 net-builder units below deliberately keep their
   minimal includes — pulling the umbrella into each would re-parse the whole engine per unit and
   regress generation time. The emitted kernel stays minimal too (runtime.hpp + sun_data.hpp). *)
        "#include \"numtracer/numtracer.hpp\"\n",
(* the two contraction phases (phase A: contract each distinct trace once over a flat work list;
   phase B: fold each net's traces as a balanced tree). Host-only (it spawns threads) and MAIN-TU
   ONLY — the net-builder units are compiled -fno-exceptions and must not see it. *)
        "#include \"numtracer/numeric/trace_fold.hpp\"\n",
        "#include <iostream>\n#include <string>\n#include <utility>\n#include <iterator>\n#include <vector>\n#include <array>\n#include <cstdlib>\n",
        "#include <thread>\n#include <atomic>\n#include <mutex>\n#include <chrono>\n#include <cstdio>\n#include <algorithm>\n#include <system_error>\n",
        "#include <unistd.h>\n#if defined(__GLIBC__)\n#include <malloc.h>\n#endif\n",
(* Resident set in MB, straight from /proc — so the generator can report its OWN peak rather than
   the caller having to wrap it in `/usr/bin/time`. Returns 0 where /proc is absent. *)
        "static double ntRssMB(){ long pages=0; if(FILE* f=std::fopen(\"/proc/self/statm\",\"r\")){ long tot=0; if(std::fscanf(f,\"%ld %ld\",&tot,&pages)!=2) pages=0; std::fclose(f); } return pages*(double)sysconf(_SC_PAGESIZE)/1048576.0; }\n",
        "#include <sstream>\n#include <unordered_map>\n",
        "using numtracer::Cx;\n",
        "namespace numtracer::network {\n",
        tmpl,
        "}\n",
        "using namespace numtracer::network;\nusing namespace numtracer::numeric;\n",
(* colour-net table: forward decls + assembler only — the chunk DEFINITIONS ride the -O0 units
   (colChunkDefs, hoisted above), keeping the 6+ MB table off the serial -O1 main-TU compile. *)
        colMainDecls];
(* The shared header block goes through a per-flow `_pch.hh` so the build can precompile it ONCE
   instead of per unit TU. Measured 2026-08-08 (za3_147, one -O0 unit, perf instructions:u):
   11.10 G plain -> 2.29 G with -include-pch, i.e. 4.85x less compile work per unit, against a
   ~2 s one-off PCH build amortised over 8 (za3_147) to 74 (za4_147) units.
   The `#ifndef NT_GEN_PCH` guard keeps the emitted source STANDALONE-compilable — required, since
   ab_gen.sh and any hand build compile these TUs with no PCH at all. The guard must suppress the
   textual include when a PCH is in play: re-including the same headers on top of the PCH costs
   5.98 G, throwing away half the win. *)
    unitInc =
      "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n#include \"numtracer/core/lit.hpp\"\n#include <utility>\n" <>
(* dressed nets emit dch<i>()/dsl<i>() builders here (the big DChainTok/DSlot literals — moved OFF
   the single -O1 main TU onto these parallel -O0 units, since a 100k+-char braced-init is ~quadratic
   even at -O0); they need the dressed-token types from numeric_contract.hpp. Non-dressed units don't
   include it (stay byte-identical + fast). *)
        If[hasDressed,
          "#include \"numtracer/numeric/numeric_contract.hpp\"\n",
          ""
        ] <>
(* colour-net chunk defs ride the units now (see colChunkDefs above); they build SUNNet literals
   through function-local SUNEnvs, so those units need the SU(N) engine header the net builders
   otherwise don't touch. Gated: colour-free flows keep their units byte-identical. *)
        If[colChunkDefs =!= {},
          "#include \"numtracer/network/sun_net.hpp\"\n",
          ""
        ];
    unitPre =
      "// GENERATED by MakeNTKernel — do not edit. Numeric net-builder unit (compiled -O0).\n" <>
        "#ifndef NT_GEN_PCH\n" <> unitInc <> "#endif\n" <>
        "using numtracer::Cx;\nnamespace numtracer::network {\n" <> tmpl <> "}\nusing namespace numtracer::network;\n" <>
        If[hasDressed,
          "using namespace numtracer::numeric;\n",
          ""];
(* NET-LEVEL CSE: dense projections (e.g. the σ^μν struct-7 quark-gluon vertex) emit the SAME net
   sub-term thousands of times — measured 48624 lnet terms but only 66 distinct (534x), so the raw
   generator C++ balloons to ~25 MB and the -O0 compile dominates generation. Hash-cons each DISTINCT
   net term into a shared accessor `lc<k>()` / `dc<k>()` (a function-local `static const` so the net
   is also BUILT once at run time, not once per use), and reference it. Trivial/empty literals and
   unique terms stay inline. Correctness-preserving: each use copies the shared NetVal/DiracNet,
   exactly as the inlined expression did. *)
    cseDefs = {};
    cseDecls = "";
(* The lc<k>()/dc<k>() net-level CSE runs for BOTH paths: the dressed (collected) lnet builders are
   just as repetitive as the dense σ case (the same Lorentz/projector structures), and the accessors
   only change the GENERATOR's internal sharing — each use copies the shared NetVal/DiracNet, so the
   contracted result (hence the committed kernel) is byte-identical. Previously skipped when dressed,
   which left the dressed generator C++ bloated and its -O0 compile ~2x slower. *)
    Module[{dPool = dsGet[], lPool = lsGet[], lCnt, dCnt, lMap = <||>, dMap = <||>, li = 0, di = 0},
(* Counted by ID, not by string: the key columns are already interned, so this is a Counts over a
   packed integer vector and a Lookup into id order. `Range[0, n-1]` recovers the counts in id order,
   which IS the first-appearance order the old `KeyValueMap` over `Counts[strings]` walked — that is
   what keeps the lc<k>/dc<k> numbering identical. *)
      With[{
        ntT =
          First @
            AbsoluteTiming[
              lCnt = Lookup[Counts[Flatten[lorentzNetStrs]], Range[0, Length[lPool] - 1], 0];
              dCnt = Lookup[Counts[Flatten[diracNetStrs]], Range[0, Length[dPool] - 1], 0];]},
        ntLog["[prof] CSE Counts (", Total[lCnt], "+", Total[dCnt], " terms): ", ntT, " s"]];
      Do[
        With[{t = lPool[[k]]},
          If[lCnt[[k]] >= 2 && t =!= "NetVal{}" && t =!= "",
            lMap[t] = "lc" <> str[li];
            li++]],
        {k, Length[lPool]}];
      Do[
        With[{t = dPool[[k]]},
          If[dCnt[[k]] >= 2 && t =!= "DiracNet{}" && t =!= "",
            dMap[t] = "dc" <> str[di];
            di++]],
        {k, Length[dPool]}];
      cseDefs =
        Join[
          KeyValueMap[
            Function[{t, nm},
              "const DiracNet& " <> nm <> "(){ static const DiracNet v = " <> t <> "; return v; }"],
            dMap],
          KeyValueMap[
            Function[{t, nm},
              "const NetVal& " <> nm <> "(){ static const NetVal v = " <> t <> "; return v; }"],
            lMap]];
      cseDecls =
        StringJoin[
          Riffle[
            Join[
              KeyValueMap[
                Function[{t, nm},
                  "const DiracNet& " <> nm <> "();"],
                dMap],
              KeyValueMap[
                Function[{t, nm},
                  "const NetVal& " <> nm <> "();"],
                lMap]],
            "\n"]];
(* The rewrite is now a POOL substitution. Every occurrence of a net term shares one pool entry, so
   replacing that entry replaces every occurrence — O(#distinct) instead of O(#sub-terms). The
   previous spelling was already vectorised (13.5 s -> 1.1 s by hoisting the substitution onto the
   distinct keys); this removes the remaining 12.7 M-element Lookup entirely. Injective, because an
   `lc<k>()` / `dc<k>()` accessor can never collide with a net literal, so the pool's first-appearance
   order — and hence every downstream id — is untouched. *)
      With[{
        ntT =
          First @
            AbsoluteTiming[
              dPoolCse = If[KeyExistsQ[dMap, #], dMap[#] <> "()", #]& /@ dPool;
              lPoolCse = If[KeyExistsQ[lMap, #], lMap[#] <> "()", #]& /@ lPool;]},
        ntLog["[prof] CSE ref-rewrite: ", ntT, " s"]];
      ntLog["[cse] net terms: lnet ", Total[lCnt], "->", Length[lMap], " distinct, dnet ", Total[dCnt], "->", Length[dMap], " distinct shared builders"]
    ];
(* ---- GLOBAL SUB-TERM DEDUP ------------------------------------------------------------------
   A net is Σ_b scal_b · contract(dn_b, ln_b, dch_b, dsl_b). The generator's cost is one trace
   contraction per (net, sub-term) — but the SAME (dn,ln,dch,dsl) tuple recurs across nets and
   colour branches, so most of those contractions recompute a trace already computed. Measured on
   the dense flows: 30,807 contractions for 6,041 distinct traces (5.1x), and 246,456 for 32,784
   (7.5x). So contract each DISTINCT trace ONCE into a shared table and let every net fold the
   table with its own scalars. Two independent wins:
     - the contraction phase (the bottleneck) shrinks by the redundancy factor;
     - the parallel phase becomes a FLAT list of uniform work items. Scheduling per NET could not
       use the machine: sub-terms per net are wildly skewed (max 2880 vs a median of 27), so the
       single biggest net alone exceeded the ideal per-thread load and pinned utilisation at ~33%
       however many cores were available.
   Sub-terms sharing a trace are merged and their scalars SUMMED (which also shortens each net's
   fold), and a merged term whose scalars sum to 0 is dropped (σ-commutator cancellations).

   This SUPERSEDES the old per-net (dn,ln) merge, which ran only when !hasDressed on the theory that
   "the collected path emits ~1 sub-term per net, so there is nothing to merge". That was false for
   dense flows (187 sub-terms/net). The global key reduces to (dn,ln) on the non-dressed path (where
   dch/dsl are constant ""), so it does everything the old merge did, and across nets as well.

   NT_GEN_NO_DEDUP=1 turns the dedup off: every occurrence becomes its own trace, nothing is merged
   or dropped, and nothing is cached (nReused=0), so each net contracts its own sub-terms on demand
   — the pre-dedup behaviour. It is the escape hatch, and the control for the equivalence test
   (generate twice, compare the kernels' VALUES; a byte-diff is meaningless because dedup changes
   GlobalEnv interning order and so renumbers every sN). *)
    ntNoDedup = ntEnvFlag["NT_GEN_NO_DEDUP"];
(* The dedup JOIN. A hash-join over EVERY sub-term — 1.6 M on ZAAqbq1, 12.7 M on ZAAqbq2 — and
   profiling (2026-08-18) put the emitter containing it at 45-54% of the whole Wolfram phase, ahead
   of the per-diagram net build and ~30x ahead of all of FunKit's COEN lowering.

   The cost was never the join; it was the KEY. Written literally, each sub-term's key is a fresh
   4-element list of net-builder STRINGS, and that list is then hashed three separate times (GatherBy
   per net, Counts over the merged terms, and the final AssociationThread lookup). So: INTERN each
   key column ONCE into integers, pack the four ids into a single integer by mixed radix, and let
   GatherBy/Counts/Lookup work on packed machine integers. Every distinct string is hashed once
   instead of millions of times, and the per-sub-term work becomes listable arithmetic.

   BYTE-IDENTITY. The emitted trace numbering must not move, and it does not: interning assigns ids
   in FIRST-APPEARANCE order, so the key SEQUENCE fed to Counts is the same permutation as before;
   `ReverseSort` on an Association is stable (verified), so equal reference counts keep that order;
   and PositionIndex groups in first-appearance order exactly as GatherBy did (verified). Confirmed
   end-to-end by regenerating the dressed fixtures before and after and diffing the kernels.

   NT_GEN_NO_DEDUP keeps its meaning: the grouping id becomes the sub-term's global occurrence index,
   so nothing ever merges. The PACKED key is kept alongside, because the emitted tables still need
   each distinct entry's four string columns; on the normal path the two are the same number.

   LEVER (b): the traceKey is dressing-free {ds, ls, chain, structural-options} — combinations that
   share a concrete structure but differ only in dressing collapse to ONE trace. The dressKey is the
   sorted dress-atom multiset; it becomes the DPoly channel a net's fold routes this sub-term's
   (scaled) trace into. Sub-terms MERGE only when they share BOTH the trace AND the dress channel (a
   merge across channels would corrupt the DPoly). For a non-dressed flow every dressKey is {}, so
   the grouping reduces to the old traceKey grouping. *)
    With[{ntT = First @ AbsoluteTiming[
      With[{joined =
          ntGenDedupJoin[
            diracNetStrs, lorentzNetStrs, subScalars, dressChains, dressSlotOpts, dressAtomIds,
            dPoolCse, lPoolCse, dcGet[], dlGet[], drGet[], hasDressed, ntNoDedup]},
        subKeysLen   = joined["subKeysLen"];
        netTraceRows  = joined["netTraceRows"];
        netDressRows  = joined["netDressRows"];
        netScalarRows = joined["netScalarRows"];
        distinctSubs = joined["distinctSubs"];
        nSub         = joined["nSub"];
        nReused      = joined["nReused"];
(* Published for the COMPILE step (mainOpt), which runs later, in another function, and needs the
   flow's distinct-trace count to pick the main-TU optimisation level — the generator RUN scales with
   it. A global rather than a threaded argument because the two are already strictly sequential
   within one generation. Set HERE, at the hand-off, so a stale value from a previous flow can never
   leak into the next one's decision. *)
        $ntGenNSub = nSub;]]},
      ntLog["[prof] sub-term dedup join: ", ntT, " s"]];
    ntLog[
      "[cse] sub-terms: ",
      Total[subKeysLen],
      " contractions -> ",
      nSub,
      " distinct traces (",
      ToString @ NumberForm[N[Total[subKeysLen] / Max[1, nSub]], {5, 2}],
      "x), ",
      nReused,
      " reused (cached), ",
      nSub - nReused,
      " singletons; per-net folds total ",
      Total[Length /@ netTraceRows],
      " (longest ",
      Max[Append[Length /@ netTraceRows, 0]],
      ")",
      If[ntNoDedup,
        " [NT_GEN_NO_DEDUP]",
        ""]];
(* The distinct traces are emitted as ONE flat table each (still chunked by ntChunkDef, whose helpers
   the bin-packer scatters across the -O0 units); the nets reference them by index. *)
(* Chunk/intern/bin-pack the net-builder tables into the -O0 unit TUs. Separated from the main()
   data tables below because the two have different levers: this half is dominated by ntChunkDef's
   dedup scans, that half by integer-to-text. *)
    With[{ntT = First @ AbsoluteTiming[
    {sdnDefs, sdnCDecl} =
      ntChunkDefs[
        "sdn",
        "std::vector<DiracNet>",
        If[nSub === 0,
          {{}},
          {distinctSubs[[All, 1]]}]];
    {slnDefs, slnCDecl} =
      ntChunkDefs[
        "sln",
        "std::vector<NetVal>",
        If[nSub === 0,
          {{}},
          {distinctSubs[[All, 2]]}]];
(* DRESSED slot tables — INTERNED (chain pool + option pool + per-sub-term index arrays).
   Expanding each structure×dressing COMBINATION into its own single-option sub-term (so phase A
   contracts them in parallel, killing the serial dress_collect) makes the chain and slot columns
   MASSIVELY redundant: a net's thousands of combinations share ONE chain (za3_147: 5 distinct over
   12101 sub-terms) and draw their DSlotOpts from a tiny per-slot option pool (35 distinct over 84672
   emissions — 99.9% redundant). Emitting the full literals per sub-term blew the generator SOURCE to
   ~9.5 MB (compile 2.7 s -> 13.5 s). Instead emit the distinct chains (`chp`) and options (`optp`)
   ONCE, and each sub-term as compact INDICES (`sdchR`: its chain index; `sdslR`: one option index per
   slot); main rebuilds sdch[k]/sdsl[k] from them in an O(nSub) loop — the exact hash-consing the
   sidx/dsc tables already use. Source ~9.5 MB -> ~0.4 MB, so the phase-A run-time win no longer costs
   a compile-time regression. Non-dressed sub-terms carry an empty chain / empty option list (the
   sdch[k].empty() numeric_value_netval fast path is preserved). *)
    {chpDefs, chpCDecl, sdchrDefs, sdchrCDecl, optpDefs, optpCDecl, sdslrDefs, sdslrCDecl} =
      If[hasDressed,
        Module[{chainStrs, combos, uChains, chainPos, uOpts, optPos, chainDefs, chainDecl, chainRefDefs, chainRefDecl, optDefs, optDecl, slotRefDefs, slotRefDecl},
          chainStrs = If[nSub === 0, {}, distinctSubs[[All, 3]]];
          combos    = If[nSub === 0, {}, distinctSubs[[All, 4]]];(* per sub-term: its option-string LIST ({} for a non-slot sub-term) *)
          uChains = DeleteDuplicates[chainStrs];
          chainPos = AssociationThread[uChains -> Range[Length[uChains]] - 1];
          uOpts = DeleteDuplicates[Flatten[combos]];
          optPos = AssociationThread[uOpts -> Range[Length[uOpts]] - 1];
          {chainDefs, chainDecl} = ntChunkDefs["chp", "std::vector<std::vector<DChainTok>>", {uChains}];
          {chainRefDefs, chainRefDecl} = ntChunkDefs["sdchR", "std::vector<int>", {ntIntStrs[chainPos /@ chainStrs]}];
          {optDefs, optDecl} = ntChunkDefs["optp", "std::vector<DSlotOpt>", {uOpts}];
          {slotRefDefs, slotRefDecl} = ntChunkDefs["sdslR", "std::vector<std::vector<int>>",
            {ntIntRow[optPos /@ #]& /@ combos}];
          {chainDefs, chainDecl, chainRefDefs, chainRefDecl, optDefs, optDecl, slotRefDefs, slotRefDecl}],
        {{}, "", {}, "", {}, "", {}, ""}];
    chunkDecls = sdnCDecl <> slnCDecl <> chpCDecl <> sdchrCDecl <> optpCDecl <> sdslrCDecl;
    allDefs = Join[cseDefs, sdnDefs, slnDefs, chpDefs, sdchrDefs, optpDefs, sdslrDefs, colChunkDefs];
(* Size-aware unit count (~$ntUnitChars per unit, 8..$ntUnitCap): the CSE accessors are many but
   small, so the old 12-defs/unit rule would emit hundreds of tiny TUs each re-parsing the shared
   decl header. Defs are packed into the units by GREEDY BIN-PACKING (largest def first, into the
   currently-smallest unit — LPT) rather than round-robin by index: def sizes are skewed, and
   round-robin left a large max/median spread that made the biggest units the makespan stragglers.
   The cap must stay above what the size target asks for, or units quietly grow back past it. Now that
   ntChunkDef bounds any SINGLE def to ~$ntDefChunk, bin-packing can actually hit the target: a single
   oversized def is no longer an irreducible floor. *)
(* Total[StringLength/@...], NOT StringLength[StringJoin[...]] — the latter materialises every def as
   one giant string just to measure it (much slower, and it allocates the lot). *)
    nUnits = Min[Min[$ntUnitCap, Max[8, Ceiling[Total[StringLength /@ allDefs] / $ntUnitChars]]], Max[1, Length[allDefs]]];
    units =
      If[allDefs === {},
        {},
        (unitPre <> StringRiffle[#, "\n"] <> "\n")& /@
          Module[{lens = StringLength /@ allDefs, order, loads = ConstantArray[0, nUnits], bin},
            order = Reverse @ Ordering[lens];
            bin = Table[With[{b = First @ Ordering[loads, 1]}, loads[[b]] += lens[[d]]; b], {d, order}];
            Lookup[GroupBy[Transpose[{bin, allDefs[[order]]}], First -> Last], Range[nUnits], {}]]];
(* the shared decl header: net-builder + CSE-accessor forward declarations (the units that call the
   lc<k>()/dc<k>() accessors #include this — emitted ONCE here, not duplicated per unit). *)
    decl =
      "// GENERATED by MakeNTKernel — do not edit. Numeric net-builder declarations.\n#pragma once\n" <> "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n#include <vector>\n" <>
        If[hasDressed,
          "#include \"numtracer/numeric/numeric_contract.hpp\"\n",
          ""
        ] <> "using namespace numtracer::network;\n" <>
        If[hasDressed,
          "using namespace numtracer::numeric;\n",
          ""
        ] <> cseDecls <> "\n" <>
(* the DISTINCT-trace tables (see the global sub-term dedup above): one flat builder each, which
   the nets index into — not one builder per net, as before the dedup. *)
        "std::vector<DiracNet> sdn0();\n" <> "std::vector<NetVal> sln0();\n" <>
        If[hasDressed,
          "std::vector<std::vector<DChainTok>> chp0();\n" <> "std::vector<int> sdchR0();\n" <>
          "std::vector<DSlotOpt> optp0();\n" <> "std::vector<std::vector<int>> sdslR0();\n",
          ""
        ] <>
(* chunk helpers of the oversized builders — the bin-packer may put a builder's helpers in a
   different unit than its assembler, so these must be declared here, not per-unit. *)
        chunkDecls;
]},
      ntLog["[prof] unit tables + bin-pack: ", ntT, " s"]];
    (* component-table init: comp[base][mu] = <MPoly builder>, skipping structural zeros. *)
    compInit =
      StringJoin @
        KeyValueMap[
          Function[{base, comps},
            StringJoin @
              MapIndexed[
                Function[{s, mu},
                  If[s === "env.zero()",
                    "",
                    "  comp[" <> str[base] <> "][" <> str[mu[[1]] - 1] <> "] = " <> s <> ";\n"]],
                comps]],
          compCpp];
(* The main TU: on a dense flow this is ~100% flat integer tables (ZAAqbq2: 99.9% of 25.7 MB), so
   this timer measures integer-to-text throughput and nothing else. *)
    With[{ntT = First @ AbsoluteTiming[
    main =
      StringJoin[
        Flatten[
          {
            "int main(int argc, char** argv){\n",
            If[ntSingleQ[], "  numtracer::codegen::emit_precision() = numtracer::codegen::EmitPrecision::Single;\n", ""],
            "  std::string decor = \"static inline\"; std::string hns = \"" <> nsInner <> "\";\n",
            "  for(int a=1;a<argc;++a){ std::string s=argv[a]; if(s==\"-d\"&&a+1<argc) decor=argv[++a]; else if(s==\"-n\"&&a+1<argc) hns=argv[++a]; }\n",
            "  const int nsym = " <> str[nsym] <> ";\n",
(* units is emitted BEFORE the LorentzEnv because the env binds both nsym and the unit groups;
   comp/atomDen and the trace entry points are then built through `env`. *)
            "  std::vector<std::vector<int>> units = {" <> StringRiffle[ntIntRow /@ unitG, ","] <> "};\n",
            "  LorentzEnv env(nsym, units);\n",
            "  std::vector<std::array<MPoly,4>> comp(" <> str[maxBase + 1] <> ", {env.zero(),env.zero(),env.zero(),env.zero()});\n",
            compInit,
            "  std::vector<std::string> symNames = {" <> StringRiffle[("\"" <> # <> "\"")& /@ symNames, ","] <> "};\n",
(* the DISTINCT-trace tables (see the global sub-term dedup). Flat, indexed by trace id: a plain
   trace has an empty chain (contract via sdn[k]), a dressed (structural) one uses sdch[k]/sdsl[k] via
   numeric_value_dressed_netval_mp. sdch/sdsl are emitted only when the kernel has dressed nets. *)
            "  std::vector<DiracNet> sdn = sdn0();\n",
            "  std::vector<NetVal> sln = sln0();\n",
            If[hasDressed,
(* rebuild the per-sub-term chain/slot tables from the interned pools (chp/optp) + index arrays
   (sdchR/sdslR) — O(nSub), reproduces the full sdch/sdsl exactly, so the trace lambda is untouched. *)
              "  std::vector<std::vector<DChainTok>> chp = chp0(); std::vector<int> sdchR = sdchR0();\n" <>
              "  std::vector<DSlotOpt> optp = optp0(); std::vector<std::vector<int>> sdslR = sdslR0();\n" <>
              "  const size_t NSD = sdchR.size();\n" <>
              "  std::vector<std::vector<DChainTok>> sdch(NSD); std::vector<std::vector<DSlot>> sdsl(NSD);\n" <>
              "  for(size_t k=0;k<NSD;++k){ sdch[k]=chp[sdchR[k]]; sdsl[k].reserve(sdslR[k].size());\n" <>
              "    for(int oi: sdslR[k]) sdsl[k].push_back(DSlot{optp[oi]}); }\n",
              ""],
(* per net: which traces it references, and with what scalar (sub-terms sharing a trace have
   already been merged, and zero sums dropped, at codegen time) *)
(* sidx/dsc HASH-CONSED. Written out in full these two tables dominate the main TU: on the
   four-quark Fierz gate they were 0.29 MB and 2.14 MB of literals for 3470 nets, and since the
   non-dressed main TU compiles at -O2 (see mainOpt) a single multi-megabyte braced-init costs
   MINUTES — measured >600 s at -O2 vs 6 s at -O0 on the same file, while every net-builder unit
   took 0.65 s. $ntDefChunk already solved this for the net-builder units; these data tables were
   never covered by it.
   The redundancy is extreme because sub-terms sharing a trace share their scalars: 73144 Cx
   literals over only 601 DISTINCT values, and 3470 index rows over 366 distinct. So dedupe rows
   for sidx, and for dsc dedupe BOTH levels (values, then the rows of value-indices — row dedup
   alone leaves ~1 MB, since the rows differ while their entries repeat). The runtime rebuild
   below is O(nets) and reproduces sidx/dsc EXACTLY as before, so fold_nets is untouched. *)
            With[{
              idxRows = netTraceRows,
              (* the per-sub-term dressing monomials, row-deduped like sidx *)
              drRows = netDressRows,
              scaRows = netScalarRows},
              With[{distinctIdxRows = DeleteDuplicates[idxRows], distinctScalars = DeleteDuplicates[Flatten[scaRows]], distinctDressMonos = DeleteDuplicates[Flatten[drRows, 1]]},
                With[{idxPos = AssociationThread[distinctIdxRows -> Range[Length[distinctIdxRows]] - 1], valPos = AssociationThread[distinctScalars -> Range[Length[distinctScalars]] - 1], drValPos = AssociationThread[distinctDressMonos -> Range[Length[distinctDressMonos]] - 1]},
                  With[{scaIdxRows = Map[valPos, scaRows, {2}], drIdxRows = Map[drValPos, drRows, {2}]},
                    With[{distinctScalarRows = DeleteDuplicates[scaIdxRows], distinctDressRows = DeleteDuplicates[drIdxRows]},
                      With[{scaPos = AssociationThread[distinctScalarRows -> Range[Length[distinctScalarRows]] - 1], drRowPos = AssociationThread[distinctDressRows -> Range[Length[distinctDressRows]] - 1]},
                        StringJoin[
                          "  const size_t NNET = " <> str[Length[netTraceRows]] <> ";\n",
(* the four big index/scalar tables move to top-level chunk functions (see ntBigTableFns): together
   with colnets they are what makes main() large enough to break the compiler. *)
                          emitBigTable["ntSidxU", "std::vector<std::vector<int>>", ntIntRow /@ distinctIdxRows,
                            "std::vector<std::vector<int>> sidxU"],
                          emitBigTable["ntSidxR", "std::vector<int>", ntIntStrs[idxPos /@ idxRows],
                            "std::vector<int> sidxR"],
                          "  std::vector<std::vector<int>> sidx(NNET);\n",
                          "  for(size_t i=0;i<NNET;++i) sidx[i]=sidxU[sidxR[i]];\n",
                          "  std::vector<Cx> dscV = {" <>
                            StringRiffle[
                              Function[s,
                                  "Cx{" <> cppNum[Re[s]] <> "," <> cppNum[Im[s]] <> "}"
                                ] /@ distinctScalars,
                              ","
                            ] <> "};\n",
                          emitBigTable["ntDscU", "std::vector<std::vector<int>>", ntIntRow /@ distinctScalarRows,
                            "std::vector<std::vector<int>> dscU"],
                          emitBigTable["ntDscR", "std::vector<int>", ntIntStrs[scaPos /@ scaIdxRows],
                            "std::vector<int> dscR"],
                          "  std::vector<std::vector<Cx>> dsc(NNET);\n",
                          "  for(size_t i=0;i<NNET;++i){ const auto& r=dscU[dscR[i]]; dsc[i].reserve(r.size());\n",
                          "    for(int k: r) dsc[i].push_back(dscV[k]); }\n",
(* LEVER (b): the per-sub-term dressing monomials sdr[i][j], deduped on BOTH levels exactly like dsc — a
   dense dressed flow has thousands of sub-terms but only a handful of DISTINCT dressing monomials (sdrV)
   and few distinct index-rows (sdrU), so a flat braced-init would be huge (measured: 115 KB single row on
   za3_147, +37 s compile) while this stays a few KB. The dressed phase-B fold routes each sub-term's scaled
   MPoly trace into its DPoly channel sdr[i][j] (empty monomial = undressed). Emitted only for dressed flows. *)
                          If[hasDressed,
                            "  std::vector<DMono> sdrV = {" <>
                              StringRiffle[ntIntRow /@ distinctDressMonos, ","] <> "};\n" <>
                            "  std::vector<std::vector<int>> sdrU = {" <>
                              StringRiffle[ntIntRow /@ distinctDressRows, ","] <> "};\n" <>
                            emitBigTable["ntSdrR", "std::vector<int>", ntIntStrs[drRowPos /@ drIdxRows],
                              "std::vector<int> sdrR"] <>
                            "  std::vector<std::vector<DMono>> sdr(NNET);\n" <>
                            "  for(size_t i=0;i<NNET;++i){ const auto& r=sdrU[sdrR[i]]; sdr[i].reserve(r.size());\n" <>
                            "    for(int k: r) sdr[i].push_back(sdrV[k]); }\n",
                            ""]]]]]]]],
(* the projector atom denominators are keyed by ATOM ID (e.inv/e.invS), not by position, and each
   id is filled idempotently — and the distinct traces cover every lnet that occurs. So scanning
   the deduped table gives the same atomDen as scanning every occurrence did. *)
            "  auto atomDen = env.collect_atom_denoms(sln, comp);\n",
            "  for(auto &a: atomDen) a = reduce_units(a, units);  // bare-loop k^2 -> monomial l1^2 -> cancels\n",
(* MATSUBARA EVENNESS, proven while contracting. If every trace and every atom denominator carries
   only EVEN powers of the Matsubara frequency, the kernel satisfies kernel(+w) == kernel(-w) and
   DiFfRG's QuadratureIntegrator_fT may collapse `kernel(+w) + kernel(-w)` to `2*kernel(w)`: half
   the Matsubara-sum work at runtime, and one fewer inlined copy of the whole kernel body per
   launch (which on a large flow is the dominant ptxas cost).

   Proven HERE, not by DiFfRG's MakeKernel "MatsubaraEven" option, which cannot work on this path:
   MakeKernel is handed the placeholder `body = 0.` (DiFfRG_compat.m) because NumTracer overwrites
   kernel.hh afterwards, so its `PossibleZeroQ[Simplify[expr - (expr /. w -> -w)]]` is trivially
   True for EVERY flow. Passing that option through would stamp the trait on kernels that are not
   even and silently drop the odd half of the sum. The numeric path has no Mathematica expression
   to test — the body is a set of polynomials this generator computes — so the proof has to live
   where the polynomials do.

   The test is a SUFFICIENT condition (odd terms that cancel between monomials read as odd), which
   is the safe direction: a false "odd" costs an optimisation, a false "even" is wrong physics. *)
            If[mIdx >= 0,
              "  // Matsubara evenness (see poly_even_in): every trace and every atom denominator must\n" <>
              "  // carry only even powers of var(" <> str[mIdx] <> "), the Matsubara frequency.\n" <>
              "  std::atomic<bool> ntMEven{true};\n" <>
              "  for(const auto &a: atomDen) if(!poly_even_in(a, " <> str[mIdx] <> ")) ntMEven.store(false, std::memory_order_relaxed);\n",
              ""],
            "  const bool ntprof = (std::getenv(\"NT_GEN_PROFILE\")!=nullptr);\n",
            "  unsigned workersA=std::thread::hardware_concurrency(); if(!workersA)workersA=4u;\n",
            "  if(const char* mw=std::getenv(\"NT_GEN_MAXW\")){int v=std::atoi(mw); if(v>0&&(unsigned)v<workersA)workersA=(unsigned)v;}\n",
(* SEPARATE worker count for phase B. The two phases have very different memory profiles per
   worker, so one knob cannot tune both. Phase A's transient is `hw` concurrent contractions, but
   it is trimmed away afterwards. Phase B's is `hw` concurrent RECOMPUTES of the uncached traces —
   and those are the SINGLETONS, which are the heavy ones (Codegen.m orders by descending refcount,
   and a trace that recurs is a simple structure while a unique one is complex: on ZAAqbq1 the 2550
   reused traces are 20.1 MB total, ~8 KB each, while the 1164 singletons could not be cached at
   all inside 10 GB). Those recomputes land on top of the live window, so phase B can need FEWER
   workers than phase A even though it is the cheaper phase in CPU terms. Defaults to `hw`. *)
            "  unsigned workersB=workersA; if(const char* mb=std::getenv(\"NT_GEN_MAXW_B\")){int v=std::atoi(mb); if(v>0)workersB=(unsigned)v;}\n",
            With[{
              (* LEVER (b): the trace table is PLAIN MPoly for BOTH paths now — a dressed sub-term's
                        structural trace is a plain MPoly (numeric_value_dressed_netval_mp) and its dressing
                        rides the per-sub-term scalar (dsc) + monomial (sdr), assembled into a DPoly only in
                        the phase-B fold. So phase A caches MPoly either way. *)
              PT = "MPoly"},
              StringJoin[
                {
                  "  const long NSUB = " <> str[nSub] <> ";\n",
(* how many traces are RESIDENT. Default: the reused ones (refCount >= 2), which the codegen-time
   ordering puts first — a singleton is contracted once whether cached or not, so caching it is
   pure RAM for no saving. NT_GEN_MEMO_MAX overrides either way (clamped to [0, NSUB]): lower it
   when memory is tight (the RAM lever — the dense flows are memory-bound before they are
   compute-bound), raise it to NSUB to put the singletons in phase A too, which costs their RAM
   but gives phase A the whole work list to balance.
   DRESSED (collected-slot) flows default to NSUB: even after lever (b) dedups the dressing variants of a
   concrete trace (so nReused is no longer ~0), the plain-MPoly traces are individually SMALL and phase B
   is parallel over NETS not traces — leaving the singletons to phase B lets the one dominant net serialise
   thousands of contractions. Caching them all (nSub) is RAM-cheap now the trace table has no dressing
   dimension, and lets phase A contract them over its flat W-parallel work list. Memory-bound dressed flows
   (ZAAqbq) dial it back with NT_GEN_MEMO_MAX. *)
                  "  long nCache = " <> str[If[hasDressed, nSub, nReused]] <> ";\n",
                  "  if(const char* mm=std::getenv(\"NT_GEN_MEMO_MAX\")){ long v=std::atol(mm); if(v>=0) nCache=std::min<long>(v,NSUB); }\n",
                  "  auto trace=[&](int k)->" <> PT <> "{\n",
(* The parity probe sits on the trace lambda rather than on the trace table T, because with
   nCache == 0 that table is EMPTY — phase B recomputes through this lambda instead. Every distinct
   trace value passes through here exactly once (a cached one was put in the cache by this same
   call), so this sees all of them and none twice. Cost is O(terms) against a contraction that is
   already superlinear in the same terms, i.e. noise.

   DRESSED flows are deliberately excluded: lever (b) strips the dressing into separate monomials
   (sdr) that the group fold multiplies back in, and those are not covered by this probe. Claiming
   evenness from the structural trace alone would be unsound if a dressing carried an odd power. *)
                  If[mIdx >= 0 && !hasDressed,
                    "    " <> PT <> " tracePoly = env.numeric_value_netval(sdn[k], sln[k], comp, atomDen);\n" <>
                    "    if(!poly_even_in(tracePoly, " <> str[mIdx] <> ")) ntMEven.store(false, std::memory_order_relaxed);\n" <>
                    "    return tracePoly;\n",
                    If[hasDressed,
                      (* structural trace → plain MPoly: a non-slot sub-term contracts via sdn[k]/sln[k],
                                   a collected one via the dressing-free _mp variant (dressing stripped at codegen). *)
                      "    return sdch[k].empty()\n" <> "      ? env.numeric_value_netval(sdn[k], sln[k], comp, atomDen)\n" <> "      : env.numeric_value_dressed_netval_mp(sdch[k], sdsl[k], sln[k], comp, atomDen);\n",
                      "    return env.numeric_value_netval(sdn[k], sln[k], comp, atomDen);\n"]],
                  "  };\n",
                  (* PHASE A — contract each distinct trace once, parallel over a FLAT work list (numeric/trace_fold.hpp). *)
                  "  auto tA=std::chrono::steady_clock::now();\n",
                  "  std::vector<" <> PT <> "> traceTable = env.contract_traces<" <> PT <> ">(nCache, workersA, trace);\n",
                  "  if(ntprof){ std::size_t tb=0; for(auto &p: traceTable) tb+=poly_bytes(p);\n",
                  "    std::fprintf(stderr,\"[num] phase A: %ld distinct traces, %ld cached, table %.1f MB, %.1f s (W=%u)\\n\",\n",
                  "      NSUB, nCache, tb/1048576.0, std::chrono::duration<double>(std::chrono::steady_clock::now()-tA).count(), workersA); }\n",
(* RELEASE PHASE A'S ARENA. Phase A contracts `hw` traces CONCURRENTLY, and a single dense 4-point
   contraction transiently allocates ~1 GB — so its working set is ~hw GB even though the table it
   leaves behind is 20 MB. glibc frees that into the per-thread arenas but does not munmap it, so
   RSS stays at the phase-A high-water mark and phase B starts from a multi-GB floor instead of
   from the table. Measured on ZAAqbq1 at W=6: RSS ~7.4 GB after phase A against a 20.1 MB table,
   which is what then pushed phase B over a 10 GB cap. malloc_trim gives it back to the OS.
   Costs milliseconds, once. *)
                  "#if defined(__GLIBC__)\n",
                  "  { const double rssPre = ntRssMB(); malloc_trim(0);\n",
                  "    if(ntprof) std::fprintf(stderr,\"[num] arena trim after phase A: RSS %.0f -> %.0f MB\\n\", rssPre, ntRssMB()); }\n",
                  "#endif\n",
(* PHASE B is NOT emitted here — it is fused into the group/lowering loop below (search
   fold_groups_streaming). It used to be `mp = env.fold_nets(...)`, one fully-expanded polynomial per
   net, ALL of them returned and then held for the rest of main() while the group loop summed them
   into a second full set of per-group accumulators. On the dense 4-point flows that is 20+ GB (488
   nets x ~41 MB) against a 20 MB trace table — the generator was OOM-killed before it could emit.
   Nothing is revisited (each net polynomial is written once, read once by its group, then dead), so
   the fix is to consume it as a stream: fold each group's nets on demand, lower the group, free it.
   Phase B therefore needs `groups`, `colv`, `g` and `realOnly`, which are only declared further
   down — hence the move. `tB` still starts here so the reported phase-B time is comparable. *)
                  "  auto tB=std::chrono::steady_clock::now();\n"}]],
(* colnets HASH-CONSED (0.57 MB -> 0.115 MB on the Fierz gate: 3470 rows, 719 distinct). Two nets
   that differ only in their Lorentz/Dirac part share a colour net, so the duplication is structural.
   This also collapses the sun_value_cx calls to the distinct nets — 4.8x fewer on that flow — which
   is a RUN saving on top of the compile one.
   The table itself is built by ntColNets() in the preamble (see there): kept inside main() it made
   the optimised TU uncompilable on a large colour graph. The SUNEnvs moved with it. *)
            With[{uCol = DeleteDuplicates[colourNets]},
              With[{colPos = AssociationThread[uCol -> Range[Length[uCol]] - 1]},
                StringJoin["  std::vector<SUNNet> colnetsU = ntColNets();\n",
                  emitBigTable["ntColR", "std::vector<int>", ntIntStrs[colPos /@ colourNets],
                    "std::vector<int> colR"], "  std::vector<numtracer::Cx> colvU(colnetsU.size());\n", "  for(size_t i=0;i<colnetsU.size();++i) colvU[i]=sun_value_cx(colnetsU[i]);\n", "  std::vector<numtracer::Cx> colv(" <> str[nNet] <> ");\n", "  for(int i=0;i<" <> str[nNet] <> ";++i) colv[i]=colvU[colR[i]];\n"]
              ]],
            emitBigTable["ntGroups", "std::vector<std::vector<int>>", ntIntRow /@ groups,
              "std::vector<std::vector<int>> groups"],
            "  // `genv`, not `env`: `env` above is the LorentzEnv (nsym + unit groups) that mints and\n" <>
            "  // contracts polynomials. This is the GLOBAL SYMBOL environment the lowering interns\n" <>
            "  // fundamental symbols into. Naming it `env` would shadow the other and silently rebind\n" <>
            "  // every env.contract_traces / env.fold_groups_streaming call below.\n" <>
            "  GlobalEnv genv;\n",
            "  std::vector<GenProg> progs;\n",
(* realOnly[gi]: this group's dressing coeff is real, so only Re(trace) is consumed -> emit a
   double trace (no dead imaginary half). Defaults to all-0 (full complex) if not supplied. *)
            "  std::vector<int> realOnly = {" <>
              StringRiffle[
                If[Length[realOnlyG] === nGrp,
                  (
                      If[TrueQ[#],
                        "1",
                        "0"]
                    )& /@ realOnlyG,
                  Table["0", {nGrp}]],
                ","
              ] <> "};\n",
(* CrossTraceCSE: accumulate every group's polynomial first, then lower them all through ONE
   shared CSE builder (to_genprog_fused) instead of one independent program per trace. Measured on
   ZAqbq1_147 Mq-in: 30,547 shared SSA instrs vs 47,558 independent (0.64x), lowering cost
   unchanged. See network::FusedProg. *)
(* PHASE B + group accumulation + lowering, FUSED into one streaming pass (numeric/trace_fold.hpp's
   fold_groups_streaming). `groups` partitions the nets, so a net's folded polynomial is needed by
   exactly one group and can die the moment that group has absorbed it. At most `gwin` group
   accumulators plus `hw` in-flight net polynomials are ever live, against every net polynomial AND
   every group accumulator before. NT_GEN_GROUP_WINDOW is the dial; gwin == nGrp reproduces the old
   residency exactly, so the pre-streaming behaviour stays reachable for A/B.

   Bit-identical to what it replaces: each group still left-folds its members in group order over the
   same fold_net results, and the sink runs on the calling thread for gi = 0,1,2,... ascending — which
   is what keeps GlobalEnv intern order and the shared CSE instruction stream unchanged. The scale is
   passed per branch rather than shared so each keeps its exact expression (poly*constant here vs
   scale_trace's constant*poly). *)
            "  long netWindow = numtracer::numeric::net_window((long)sidx.size(), workersB);\n",
(* UNCONDITIONAL, deliberately. `groups` must PARTITION the nets: a duplicate means a net is folded
   into the kernel twice, a gap means one is silently dropped. Either is a wrong kernel with no other
   symptom. This used to be emitted behind `if(ntprof)` — i.e. it ran only under NT_GEN_PROFILE,
   which no production regeneration sets, so the invariant has effectively never been checked. It is
   O(nNet) once per generation (nets are tens to low thousands) against a kernel build measured in
   seconds to minutes, so gating it on a profiling flag bought nothing and cost the guard entirely.
   It is FATAL (std::exit(1) — see check_group_partition in numeric/trace_fold.hpp), and runs on
   every generation. It was made fatal after a full regeneration of all 29 DEFAULT_FLOWS reported
   zero violations, so no committed flow legitimately produces a non-partition. *)
            "  numtracer::numeric::check_group_partition(groups, " <> str[nNet] <> ");\n",
(* LEVER (b): the dressed fold reads a PLAIN-MPoly trace table T + the per-sub-term dressing monomials sdr,
   building the per-net DPoly channel-by-channel (fold_groups_streaming_dressed). The colour scale (scaleCx)
   and the sink are unchanged — only the per-net fold's trace type differs, so the group DPoly the sink
   receives is value-identical to the pre-lever-(b) DPoly-trace-table path. *)
            Which[
              crossCSE && hasDressed,
                "  FusedStream fstream(genv, realOnly);\n" <> "  env.fold_groups_streaming_dressed(sidx, dsc, sdr, groups, traceTable, nCache, workersB, netWindow, trace,\n" <> "    [&](int d, DPoly &&m){ return scaleCx(m, colv[d]); },\n" <> "    [&](size_t, DPoly &&acc){ fstream.add(acc); });\n" <> "  std::vector<FusedProg> fused = fstream.finish();\n",
              crossCSE,
                "  FusedStream fstream(genv, realOnly);\n" <> "  env.fold_groups_streaming<MPoly>(sidx, dsc, groups, traceTable, nCache, workersB, netWindow, trace,\n" <> "    [&](int d, MPoly &&m){ return m*env.constant(colv[d]); },\n" <> "    [&](size_t, MPoly &&acc){ fstream.add(acc); });\n" <> "  std::vector<FusedProg> fused = fstream.finish();\n",
              hasDressed,
                "  env.fold_groups_streaming_dressed(sidx, dsc, sdr, groups, traceTable, nCache, workersB, netWindow, trace,\n" <> "    [&](int d, DPoly &&m){ return scaleCx(m, colv[d]); },\n" <> "    [&](size_t gi, DPoly &&acc){ progs.push_back(to_genprog(acc, genv, realOnly[gi]!=0)); });\n",
              True,
                "  env.fold_groups_streaming<MPoly>(sidx, dsc, groups, traceTable, nCache, workersB, netWindow, trace,\n" <> "    [&](int d, MPoly &&m){ return m*env.constant(colv[d]); },\n" <> "    [&](size_t gi, MPoly &&acc){ progs.push_back(to_genprog(acc, genv, realOnly[gi]!=0)); });\n"
            ],
            "  if(ntprof) std::fprintf(stderr,\"[num] phase B+lower: %d nets in %d groups, window %ld, %.1f s (W=%u)\\n\", " <> str[nNet] <> ", " <> str[nGrp] <> ", netWindow, std::chrono::duration<double>(std::chrono::steady_clock::now()-tB).count(), workersB);\n",
(* the trace table is dead once every group has folded; emission below only needs the lowered
   instruction streams, which are orders of magnitude smaller. *)
            (* LEVER (b): the trace table T is plain MPoly for BOTH paths now. *)
            "  { std::vector<MPoly> dead; traceTable.swap(dead); }\n",
(* Emission was the one untimed stage of the run: it renders, dedups and writes every trace body
   serially, which on a multi-MB header is minutes, not noise. Time it so the [num] trail covers
   the whole run. *)
            "  const auto tEmit = std::chrono::steady_clock::now();\n",
            "  FillFormulas fm;\n",
            "  fm.var = [](int id)->std::string{\n",
            Table["    if(id==" <> str[i - 1] <> ") return \"" <> varFill[[i]] <> "\";\n", {i, 1, Length[varFill]}],
            "    return \"" <> ntZeroLit[] <> "\"; };\n",
            "  fm.inv = [&](int id)->std::string{ return \"" <> If[ntSingleQ[], "1.f", "1.0"] <> "/(\" + mpoly_to_cpp(atomDen[(size_t)id], symNames) + \")\"; };\n",
(* dressing fill (kind-2 `dress` leaves): the kernel BODY evaluates each dressing atom (where the
   regulators REG::* and the interpolator parameters are in scope) into `dr_<id>` and passes the
   VALUE to fill(); fill just stores it. So fm.dress returns the passed-in argument name. *)
            If[hasDressed,
              "  fm.dress = [](int id)->std::string{ return \"dr_\" + std::to_string(id); };\n",
              ""],
            "  std::cout << \"// GENERATED by gen_" <> nsInner <> ".cpp — do not edit.\\n\";\n",
(* Complex trace values carry an OVERRIDABLE type: nvcc lowers std::complex arithmetic through
   gcc _Complex builtins that device code silently computes as 0 (re-confirmed 2026-08-08 on
   CUDA 12.9: the whole kernel returned exactly 0.0 on device, correct values on host). A device
   consumer #defines NT_TRACE_COMPLEX (e.g. to cuda::std::complex<double>) BEFORE including the
   traces header; host consumers get std::complex<double> unchanged. The alias lives in the
   per-kernel namespace so multiple traces headers coexist in one TU. *)
            "  std::cout << \"#pragma once\\n#include <cmath>\\n" <>
              If[complexQ,
                "#include <complex>\\n",
                ""
              ] <> "namespace " <> kns <> " { namespace \" << hns << \" {\\n\";\n",
            If[complexQ,
              "  std::cout << \"#ifndef NT_TRACE_COMPLEX\\n#define NT_TRACE_COMPLEX std::complex<" <> $ntRealT <> ">\\n#endif\\nusing nt_complex_t = NT_TRACE_COMPLEX;\\n\";\n",
              ""],
(* Unqualified fma/sqrt in a float header would otherwise bind the global double ::fma/::sqrt on
   the host, silently running that arithmetic in double. *)
            If[ntSingleQ[], "  std::cout << \"using std::fma;\\nusing std::sqrt;\\n\";\n", ""],
            "  std::cout << \"template<int N> \" << decor << \" " <> $ntRealT <> " powr(" <> $ntRealT <> " x){ " <> $ntRealT <> " r=" <> If[ntSingleQ[], "1.f", "1.0"] <> "; for(int i=0;i<N;++i) r*=x; return r; }\\n\";\n",
            "  emit_env_layout(std::cout, genv);\n",
            "  std::cout << \"static inline constexpr int nenv = \" << genv.syms.size() << \";\\n\";\n",
(* The proven verdict, as a compile-time constant the kernel class picks up. Always emitted for a
   finite-T flow (rather than only when true) so the header records what was checked: a `false`
   here is a positive statement that the traces were tested and are not even, not a gap. *)
            If[mIdx >= 0,
              "  std::cout << \"// Matsubara evenness of the traces in var(" <> str[mIdx] <> "), proven from the monomial\\n\";\n" <>
              "  std::cout << \"// exponents at generation time (see poly_even_in). Consumed by the kernel class as\\n\";\n" <>
              "  std::cout << \"// DiFfRG's `matsubara_even` trait, which halves the Matsubara-sum evaluations.\\n\";\n" <>
              "  std::cout << \"static inline constexpr bool matsubara_even = \" << (" <>
                If[hasDressed, "false", "ntMEven.load(std::memory_order_relaxed)"] <>
                " ? \"true\" : \"false\") << \";\\n\";\n" <>
              "  if(ntprof) std::fprintf(stderr,\"[num] matsubara_even = %s\\n\", " <>
                If[hasDressed, "\"false (dressed flow: not checked)\"", "(ntMEven.load(std::memory_order_relaxed)?\"true\":\"false\")"] <> ");\n",
              ""],
            "  emit_fill(std::cout, genv, \"fill\", \"" <> fillArgSig <> "\", fm, decor);\n",
(* TRACE-BODY DEDUP: the grouping key is the dressing COEFFICIENT (diagData), finer than the trace
   STRUCTURE — so flows with many Feynman graphs that share a kinematic trace but differ only in
   their dressing/coupling coefficient (e.g. the σ-Yukawa hSigL: 503 groups, only 157 distinct trace
   bodies) emit hundreds of byte-identical trN. Render each trace, key on the name-independent body,
   and emit a duplicate as a one-line forwarder `trN(f){ return trK(f); }` (canonical K = first with
   that body). The forwarder's RETURN TYPE is read back off the emitted signature (the token right
   before " tr<i>(") rather than sniffed for a spelling: emit_cpp writes the complex return type as
   the `nt_complex_t` alias, so the old `find("std::complex<double> tr")` test silently never matched
   and every complex duplicate got a `double` forwarder — a hard compile error in any flow that has
   both complex traces and duplicate bodies. Only the type is taken from `s`; the decorator stays
   `decor`, so a canonical body that eff_decor out-of-lined does not drag `noinline` onto the
   one-line forwarder. The body is computed identically once per call site regardless (GCC CSEs the inlined
   identical traces — see ZA4 fusion notes), so this is a pure SOURCE-SIZE/compile win, runtime-neutral.
   Flows with no shared trace structure (ZAqbq1/4/7_147: 108/108 distinct) never hit `seen` ⇒ the
   emitted bytes are unchanged. *)
(* fused: ONE trace_all(f, t[]) instead of nGrp trN(). The body-dedup below is meaningless then
   (there is a single body), and the kernel reads tarr[i] rather than calling tr_i. *)
            If[crossCSE,
              "  emit_cpp_fused(std::cout, fused, \"trace_all\", decor);\n",
              "  { std::unordered_map<std::string,std::string> seen; seen.reserve((size_t)" <> str[nGrp] <> ");\n" <> "    for(int i=0;i<" <> str[nGrp] <> ";++i){\n" <> "      const std::string nm = \"tr\"+std::to_string(i);\n" <> "      std::ostringstream os; emit_cpp(os, progs[i], nm, decor);\n" <> "      std::string s = os.str(); std::string body = s.substr(s.find('{'));\n" <> "      auto it = seen.find(body);\n" <> "      if(it==seen.end()){ seen.emplace(std::move(body), nm); std::cout << s; }\n" <> "      else { const std::string sig = s.substr(0, s.find(\" \"+nm+\"(\")); const std::string rt = sig.substr(sig.rfind(' ')+1);\n" <> "        std::cout << decor << \" \" << rt << \" \" << nm << \"(const " <> $ntRealT <> " *f) { return \" << it->second << \"(f); }\\n\"; } } }\n"
            ],
            "  std::cout << \"}} // namespace " <> kns <> "::\" << hns << \"\\n\";\n",
            "  if(ntprof) std::fprintf(stderr,\"[num] emission: %.1f s\\n\", std::chrono::duration<double>(std::chrono::steady_clock::now()-tEmit).count());\n",
            "  return 0;\n}\n"}]];
(* The table builders are collected DURING the StringJoin above (each emitBigTable stuffs the bag as
   its argument is evaluated), so they can only be prepended here — referencing the bag INSIDE that
   StringJoin would splice in whatever was collected so far, which at that point is nothing. *)
    main = StringJoin[Internal`BagPart[tableBag, All]] <> main;]},
      ntLog["[prof] main() data tables: ", ntT, " s"]];
    {pre, units, decl, main}];
