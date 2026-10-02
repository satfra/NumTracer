(* CodegenGenerator.m: emitNumericGenerator and its sub-term dedup join (ntGenDedupJoin). Emits the
   C++ generator program that builds each net's Dirac/Lorentz/colour structure, contracts the traces
   numerically, and PRINTS the committed straight-line kernel header.
   Loaded by NumTracer.m via ntLoadPart, in the NumTracer`Private` context. *)

(* ---- STAGE: the global sub-term dedup JOIN ----------------------------------------------------
   Input: the six already-interned id columns (one ragged list per net) and the five pools they index.
   Output: the distinct traces, the per-net fold entries that reference them, and the counts the
   caching decision needs.

   WHY. A net is Σ_b scal_b · contract(trace_b), and the same trace recurs across nets and colour
   branches (5-7x on the dense flows). So each DISTINCT trace is contracted once into a shared table
   that every net folds with its own scalars. This also turns the parallel contraction into a flat
   list of uniform work items; per-net scheduling stalls on a few huge nets.

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
         packed key, so nothing merges. The substitution and the traceKeyFlat read-back are ONE
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
(* The key columns arrive already interned (ntGenExpandSubTerms), so nothing is hashed
   here: the join is integer arithmetic over ids. *)
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

(* Per net: merge the sub-terms that share a trace AND a dress channel (merging across channels would
   corrupt the DPoly), summing scalars; drop the zero sums. Do this BEFORE counting references: a
   repeat inside ONE net is a single reference, and a cancelled channel must not be contracted. Each
   net yields three columns {trace key, dress id, summed scalar} (ntMergeNetTerms). *)
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
    (* #refCount >= 2 == the length of the sorted prefix worth caching *)
    nReused = Total[UnitStep[Values[refCount] - 2]];
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
   they share. Returns {pre, units, decl, main, nSub} (nSub = distinct sub-terms, for the compile
   log). The program it emits is what actually contracts the traces and PRINTS the committed
   straight-line kernel header; nothing here contracts anything.

   emitNumericGenerator is the driver; the stages, in order:

     1. ntGenExpandSubTerms — SUB-TERM EXPANSION + INTERNING. Walk every net's cores, expand the
        dressed slot options into their Cartesian product, and intern the five key columns (Dirac net,
        Lorentz net, dressed chain, slot-option tuple, dressing-atom multiset) through ntMkIntern.
        Output: six ragged integer columns, one entry per (net, sub-term), plus the five pools.
     2. ntGenColourTables — the distinct colour nets: chunk builders for the unit TUs, the main-TU
        assembler, and the per-net index colR into them.
     3. ntGenPreambles — the main-TU and unit-TU headers.
     4. ntGenNetCSE — net-builder strings that recur across sub-terms become shared `lc<k>()` /
        `dc<k>()` accessors, and the pools are rewritten to call them.
     5. ntGenDedupJoin (above) — merge sub-terms sharing a trace AND a dress channel, count
        references, order the distinct traces by descending reference count. Output: the
        distinct-trace table and the per-net fold entries that index it.
     6. ntGenUnitSources — chunk the distinct-trace tables (ntGenTraceTables, ntGenDressedSlotTables),
        LPT bin-pack every definition into the -O0 unit TUs (ntLptBinPack), emit the decl header.
     7. main(): ntGenMainPrologue, ntGenMainTables (per-net tables, hash-consed by ntHashConsRows),
        ntGenMainPhaseA (contract each distinct trace), ntGenMainPhaseB (fold nets, drain groups,
        lower to SSA), ntGenMainEmission (print the kernel header). The big tables are top-level
        chunk functions (ntGenBigTable) whose definitions the driver prepends to main().

   WHAT MAKES THE OUTPUT REPRODUCIBLE. Every stage boundary above fixes an ORDER, and the emitted
   kernel's slot numbering follows it, so a reordering here yields a kernel that is still correct and
   still passes every oracle while differing byte-for-byte from the committed one. The four specific
   invariants are stated at ntGenDedupJoin (I1-I4). Two more fix the generator source: in stage 6
   the definitions are Joined in a fixed order and `Ordering` breaks ties by index, so that Join
   order decides which definition lands in which unit TU; in stage 7 the driver's StringJoin fixes
   the order of the table functions in front of main().

   The `Tuples` in stage 1 runs LAST-SLOT-FASTEST, and stages 1, 5 and 7 all rely on that alignment. *)

(* ---- helpers --------------------------------------------------------------------------------- *)

(* rows -> distinct rows (first-appearance order) and each row's 0-based index into them. With
   "Values", first intern the row ENTRIES the same way, then hash-cons the rows of value indices:
   rows alone differ while their entries repeat.
   Distinctness must be the Association's (exact) key equality, not SameQ: SameQ treats machine numbers
   one ulp apart as equal (0.5 vs 0.5000000000000001), so DeleteDuplicates would drop one and its
   lookup would come back Missing[KeyAbsent, ...] into the C++. Keys@PositionIndex is exact and keeps
   first-appearance order. *)
ntHashConsRows[rows_List] :=
  With[{u = Keys[PositionIndex[rows]]},
    <|"Rows" -> u, "RowIdx" -> AssociationThread[u -> Range[Length[u]] - 1] /@ rows|>];

ntHashConsRows[rows_List, "Values"] :=
  With[{vals = Keys[PositionIndex[Flatten[rows, 1]]]},
    Append[ntHashConsRows[Map[AssociationThread[vals -> Range[Length[vals]] - 1], rows, {2}]],
      "Values" -> vals]];

(* One big main() table as top-level chunk functions (ntBigTableFns): inline, these tables make
   main() large enough to break the compiler. Returns {definitions, the `<lhs> = <call>;` line}. *)
ntGenBigTable[nm_String, ret_String, rows_List, lhs_String] :=
  With[{t = ntBigTableFns[nm, ret, rows]},
    {t[[1]], "  " <> lhs <> " = " <> t[[2]] <> ";\n"}];

(* LPT bin-packing: largest definition first, into the currently-smallest unit. Def sizes are
   skewed, and ntChunkDef bounds any single def to ~$ntDefChunk, so the size target is reachable.
   `Ordering` breaks ties by index, so the input order decides the packing. Returns nUnits lists. *)
ntLptBinPack[defs_List, nUnits_Integer] :=
  Module[{lens = StringLength /@ defs, order, loads = ConstantArray[0, nUnits], bin},
    order = Reverse @ Ordering[lens];
    bin = Table[With[{b = First @ Ordering[loads, 1]}, loads[[b]] += lens[[d]]; b], {d, order}];
    Lookup[GroupBy[Transpose[{bin, defs[[order]]}], First -> Last], Range[nUnits], {}]];

(* ---- STAGE 1: sub-term expansion + interning --------------------------------------------------
   Per net, a colour group is a SUM of sub-terms: coreNets[i] = {core_b…} (a DiracNet literal for a
   gamma branch, a Lorentz NetVal for a gamma-free one), restScalars[i] = {{rest_b, scal_b}…} parallel.
   Each branch yields {ds, ls, scal, dc, dl, dr} columns (usually one sub-term). A DRESSED branch
   (ntDressedCore[chainStr, slotsStr]) expands the Cartesian product of its chain's slot options,
   each a triple {structStr, num, dr}: the structures (dl) form the trace key, ∏ num folds into the
   scalar, and the dress-atom multiset (dr) becomes the DPoly key. Non-dressed branches carry an
   empty dress key. Traces are dressing-stripped, so combinations differing only in dressing share
   ONE trace; the dressing rides the per-sub-term scalar + monomial into the phase-B fold.

   The five key columns are interned HERE, once per distinct value, so the dedup join works on
   integers. Byte-identity: `ntMkIntern` assigns ids in first-appearance order along the walk
   nets -> cores -> combinations (`Tuples` last-slot-fastest), which is the order the join relies
   on (I1). *)
ntGenExpandSubTerms[coreNets_, restScalars_] :=
  Module[{internDiracNet, diracNetPool, internLorNet, lorNetPool, internChain, chainPool,
          internSlotTuple, slotTuplePool, internDressMono, dressMonoPool, slotComboCache = <||>,
          slotCombo, cols, ntT},
    {internDiracNet, diracNetPool} = ntMkIntern[];
    {internLorNet, lorNetPool} = ntMkIntern[];
    {internChain, chainPool} = ntMkIntern[];
    {internSlotTuple, slotTuplePool} = ntMkIntern[];
    {internDressMono, dressMonoPool} = ntMkIntern[];
(* The Cartesian expansion of one core's slot options depends on `slotOpts` ALONE, so it is memoised
   per distinct option set. Returns {dlIds, drIds, nums, n} from one call so the columns can never
   desync: `Tuples` and `Flatten[Outer[...]]` both vary the last slot fastest, and that alignment
   between a structure tuple and its numeric coefficient is load-bearing. *)
    slotCombo[slotOpts_] :=
      Lookup[slotComboCache, Key[slotOpts],
        slotComboCache[slotOpts] =
          Module[{structs, nums, dress, cs, n},
            structs = #[[All, 1]]& /@ slotOpts;
            nums    = #[[All, 2]]& /@ slotOpts;
            dress   = #[[All, 3]]& /@ slotOpts;
            cs = Tuples[structs];
            n  = Length[cs];
            {Developer`ToPackedArray[internSlotTuple /@ cs],
(* by far the common case: no slot option on this core carries a dressing atom, so every
   combination's union is empty. Checked once per distinct option set. *)
             If[AllTrue[dress, AllTrue[#, # === {}&]&],
               ConstantArray[internDressMono[{}], n],
               Developer`ToPackedArray[internDressMono /@ (Sort[Catenate[#]]& /@ Tuples[dress])]],
             Flatten[Outer[Times, Sequence @@ nums]],
             n}]];
    ntT = First @ AbsoluteTiming[
      cols =
        Transpose @
          MapThread[
            Function[{cores, rss},
              If[cores === {},
                {{}, {}, {}, {}, {}, {}},
(* COLUMN-ORIENTED: a dressed flow has millions of sub-terms, so nothing is evaluated per sub-term
   at Mathematica level. Every column is a constant repeated n times, a `Tuples`, or an `Outer`
   product (see slotCombo); `Join @@@ Transpose[...]` regroups the per-core column tuples into six
   whole-net columns. *)
                Join @@@ Transpose @
                    MapThread[
                      Function[{nv, rv, scal},
                        Module[{lsStr = If[rv === "", "NetVal{}", rv]},
                          Which[
                            (* dressed numerator: expand slot options → structural sub-terms *)
                            MatchQ[nv, _ntDressedCore],
                              With[{chain = nv[[1]], slotOpts = nv[[2]]},
                                If[slotOpts === {},
                                  {{internDiracNet["DiracNet{}"]}, {internLorNet[lsStr]}, ntPackCx[{scal}], {internChain[chain]}, {internSlotTuple[{}]}, {internDressMono[{}]}},
(* dl = this combination's STRUCTURAL option-string LIST (one dressing-free structStr per chain
   slot), interned so the table emitter still pools the distinct structures; the numeric Cx folds
   into the sub-term scalar and the dress ids become the DPoly key. *)
                                  Module[{tb = slotCombo[slotOpts], n},
                                    n = tb[[4]];
                                    {ConstantArray[internDiracNet["DiracNet{}"], n],
                                     ConstantArray[internLorNet[lsStr], n],
                                     ntPackCx[scal * tb[[3]]],
                                     ConstantArray[internChain[chain], n],
                                     tb[[1]],
                                     tb[[2]]}]]],
                            (* gamma branch: DiracNet + projector rest *)
                            StringStartsQ[nv, "DiracNet"],
                              {{internDiracNet[nv]}, {internLorNet[lsStr]}, ntPackCx[{scal}], {internChain["std::vector<DChainTok>{}"]}, {internSlotTuple[{}]}, {internDressMono[{}]}},
                            (* gamma-free branch: whole net is the rest *)
                            True,
                              {{internDiracNet["DiracNet{}"]}, {internLorNet[nv]}, ntPackCx[{scal}], {internChain["std::vector<DChainTok>{}"]}, {internSlotTuple[{}]}, {internDressMono[{}]}}]]],
                      {cores, rss[[All, 1]], rss[[All, 2]]}]]],
            {coreNets, restScalars}];];
    ntLog["[prof] sub-term expansion: ", ntT, " s"];
    ntStageResult["ntGenExpandSubTerms",
      {"DiracNetIds", "LorNetIds", "SubScalars", "ChainIds", "SlotTupleIds", "DressMonoIds",
       "DiracNetPool", "LorNetPool", "ChainPool", "SlotTuplePool", "DressMonoPool"},
      <|"DiracNetIds" -> cols[[1]], "LorNetIds" -> cols[[2]], "SubScalars" -> cols[[3]],
        "ChainIds" -> cols[[4]], "SlotTupleIds" -> cols[[5]], "DressMonoIds" -> cols[[6]],
        "DiracNetPool" -> diracNetPool[], "LorNetPool" -> lorNetPool[], "ChainPool" -> chainPool[],
        "SlotTuplePool" -> slotTuplePool[], "DressMonoPool" -> dressMonoPool[]|>]];

(* ---- STAGE 2: colour-net tables ---------------------------------------------------------------
   The distinct colour nets are one `SUNNet{...}` literal each; on a large colour graph the table can
   dominate the main TU, which is the one serial -O1 compile. So the chunk functions ride the
   LPT-packed -O0 units like every other builder (ChunkDefs), with external (non-static) linkage,
   and only forward decls + the assembler ntColNets() stay in the main TU (MainDecls).
   colnets are HASH-CONSED: nets differing only in their Lorentz/Dirac part share a colour net, so
   sun_value_cx runs once per DISTINCT net and colR maps each net to it (ColRDef, MainText). *)
ntGenColourTables[colourNets_, nNet_] :=
  Module[{hc = ntHashConsRows[colourNets], uCol, mainDecls, chunkDefs, colR},
    uCol = hc["Rows"];
    {mainDecls, chunkDefs} =
      If[uCol === {},
        {"static std::vector<SUNNet> ntColNets(){ return {}; }\n", {}},
        Module[{envDecl, chunks},
          envDecl = StringJoin["SUNEnv sun" <> # <> "(" <> # <> "); "& /@ DeleteDuplicates @ Flatten @ StringCases[uCol, "sun" ~~ r : DigitCharacter.. ~~ "." :> r]];
          chunks = ntSplitByChars[uCol, 2];
          {StringJoin[
             MapIndexed["void ntColNets_c" <> ToString[#2[[1]] - 1] <> "(std::vector<SUNNet>& o);\n"&, chunks],
             "static std::vector<SUNNet> ntColNets(){ std::vector<SUNNet> o; o.reserve(" <> ToString[Length[uCol]] <> "); " <> StringJoin[Table["ntColNets_c" <> ToString[k - 1] <> "(o); ", {k, Length[chunks]}]] <> "return o; }\n"],
           MapIndexed["void ntColNets_c" <> ToString[#2[[1]] - 1] <> "(std::vector<SUNNet>& o){ " <> envDecl <> StringJoin[("o.push_back(" <> # <> "); ")& /@ #1] <> "}\n"&, chunks]}]];
    colR = ntGenBigTable["ntColR", "std::vector<int>", ntIntStrs[hc["RowIdx"]], "std::vector<int> colR"];
    ntStageResult["ntGenColourTables", {"MainDecls", "ChunkDefs", "ColRDef", "MainText"},
      <|"MainDecls" -> mainDecls, "ChunkDefs" -> chunkDefs, "ColRDef" -> colR[[1]],
        "MainText" ->
          StringJoin["  std::vector<SUNNet> colnetsU = ntColNets();\n", colR[[2]],
            "  std::vector<numtracer::Cx> colvU(colnetsU.size());\n",
            "  for(size_t i=0;i<colnetsU.size();++i) colvU[i]=sun_value_cx(colnetsU[i]);\n",
            "  std::vector<numtracer::Cx> colv(" <> ToString[nNet] <> ");\n",
            "  for(int i=0;i<" <> ToString[nNet] <> ";++i) colv[i]=colvU[colR[i]];\n"]|>]];

(* ---- STAGE 3: preambles -----------------------------------------------------------------------
   Pre: the main TU's includes, helpers and the colour-net assembler. UnitPre: the header every -O0
   net-builder unit starts with. *)
ntGenPreambles[hasDressed_, colMainDecls_String, colourUnitsQ_] :=
  Module[{tmpl, unitInc},
    (* shared wrapper templates so the net-builder strings compile. *)
    tmpl = "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal tproj(){ return projT(Mu,Nu,Lb,Inv); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal lproj(){ return projL(Mu,Nu,Lb,Inv); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int InvS> NetVal mproj(){ return projM(Mu,Nu,Lb,InvS); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int Inv,int InvS> NetVal eproj(){ return projE(Mu,Nu,Lb,Inv,InvS); }\n" <> "template<int Mu,int Nu> NetVal lmetric(){ return met(Mu,Nu); }\n" <> "template<int Lbl,int Base,int Mask> NetVal lvec(){ return vec(Lbl,Base); }\n" <> "template<int A,int B,int C,int D> NetVal leps(){ return epsilon(A,B,C,D); }\n" <> "inline NetVal konst(double c){ return NetVal{PTerm{Cx{c,0}, {}}}; }\n" <> "template<class L> struct litco;\n" <> "template<numtracer::Cx C> struct litco<numtracer::Lit<C>>{ static constexpr numtracer::Cx v=C; };\n" <> "template<class L> NetVal sc(NetVal x){ return scale(litco<L>::v, std::move(x)); }\n";
(* The shared header block goes through a per-flow `_pch.hh` so the build can precompile it ONCE
   instead of per unit TU (unit compiles are header-bound; the PCH cuts their work ~5x).
   The `#ifndef NT_GEN_PCH` guard keeps the emitted source STANDALONE-compilable (ab_gen.sh and hand
   builds use no PCH), and must suppress the textual include when a PCH is in play: re-including the
   headers on top of it throws away half the win. *)
    unitInc =
      "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n#include \"numtracer/core/lit.hpp\"\n#include <utility>\n" <>
(* dressed nets put their big DChainTok/DSlot literal builders on these parallel -O0 units (a huge
   braced-init in the serial main TU compiles ~quadratically); they need the dressed-token types from
   numeric_contract.hpp. Non-dressed units don't include it. *)
        If[hasDressed,
          "#include \"numtracer/numeric/numeric_contract.hpp\"\n",
          ""
        ] <>
(* colour-net chunk defs ride the units (stage 2); they build SUNNet literals through function-local
   SUNEnvs, so those units need the SU(N) engine header the net builders otherwise don't touch.
   Gated: colour-free flows keep their units byte-identical. *)
        If[colourUnitsQ,
          "#include \"numtracer/network/sun_net.hpp\"\n",
          ""
        ];
    ntStageResult["ntGenPreambles", {"Pre", "UnitPre"},
      <|"Pre" ->
          StringJoin[
            "// GENERATED by MakeNTKernel — do not edit. Numeric matrix-product tensor traces.\n",
(* one umbrella header pulls the whole engine API (network/dirac/gen/sun_net/numeric_*/core) for
   this single generator TU. The parallel -O0 net-builder units deliberately keep their minimal
   includes — pulling the umbrella into each would re-parse the whole engine per unit and regress
   generation time. The emitted kernel stays minimal too (runtime.hpp + sun_data.hpp). *)
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
            colMainDecls],
        "UnitPre" ->
          "// GENERATED by MakeNTKernel — do not edit. Numeric net-builder unit (compiled -O0).\n" <>
            "#ifndef NT_GEN_PCH\n" <> unitInc <> "#endif\n" <>
            "using numtracer::Cx;\nnamespace numtracer::network {\n" <> tmpl <> "}\nusing namespace numtracer::network;\n" <>
            If[hasDressed,
              "using namespace numtracer::numeric;\n",
              ""]|>]];

(* ---- STAGE 4: net-level CSE -------------------------------------------------------------------
   Dense projections emit the SAME net sub-term thousands of times (e.g. the σ^μν quark-gluon vertex:
   ~500x), ballooning the generator source and its -O0 compile. Hash-cons each recurring net term
   into a shared accessor `lc<k>()` / `dc<k>()` (a function-local `static const`, so it is also BUILT
   once at run time). Trivial/empty literals and unique terms stay inline. Each use copies the shared
   NetVal/DiracNet, so the committed kernel is unchanged. *)
ntGenNetCSE[dPool_List, lPool_List, diracNetIds_, lorNetIds_] :=
  Module[{lCnt, dCnt, lMap = <||>, dMap = <||>, li = 0, di = 0, dPoolCse, lPoolCse, ntT},
(* Counted by interned ID. `Range[0, n-1]` yields the counts in id order, i.e. first-appearance
   order, which fixes the lc<k>/dc<k> numbering. *)
    ntT = First @ AbsoluteTiming[
      lCnt = Lookup[Counts[Flatten[lorNetIds]], Range[0, Length[lPool] - 1], 0];
      dCnt = Lookup[Counts[Flatten[diracNetIds]], Range[0, Length[dPool] - 1], 0];];
    ntLog["[prof] CSE Counts (", Total[lCnt], "+", Total[dCnt], " terms): ", ntT, " s"];
    Do[
      With[{t = lPool[[k]]},
        If[lCnt[[k]] >= 2 && t =!= "NetVal{}" && t =!= "",
          lMap[t] = "lc" <> ToString[li];
          li++]],
      {k, Length[lPool]}];
    Do[
      With[{t = dPool[[k]]},
        If[dCnt[[k]] >= 2 && t =!= "DiracNet{}" && t =!= "",
          dMap[t] = "dc" <> ToString[di];
          di++]],
      {k, Length[dPool]}];
(* The rewrite is a POOL substitution: every occurrence of a net term shares one pool entry, so this
   is O(#distinct). Injective, because an `lc<k>()` / `dc<k>()` accessor can never collide with a net
   literal, so the pool's first-appearance order — and every downstream id — is untouched. *)
    ntT = First @ AbsoluteTiming[
      dPoolCse = If[KeyExistsQ[dMap, #], dMap[#] <> "()", #]& /@ dPool;
      lPoolCse = If[KeyExistsQ[lMap, #], lMap[#] <> "()", #]& /@ lPool;];
    ntLog["[prof] CSE ref-rewrite: ", ntT, " s"];
    ntLog["[cse] net terms: lnet ", Total[lCnt], "->", Length[lMap], " distinct, dnet ", Total[dCnt], "->", Length[dMap], " distinct shared builders"];
    ntStageResult["ntGenNetCSE", {"Defs", "Decls", "DiracPool", "LorPool"},
      <|"Defs" ->
          Join[
            KeyValueMap["const DiracNet& " <> #2 <> "(){ static const DiracNet v = " <> #1 <> "; return v; }"&, dMap],
            KeyValueMap["const NetVal& " <> #2 <> "(){ static const NetVal v = " <> #1 <> "; return v; }"&, lMap]],
        "Decls" ->
          StringJoin[
            Riffle[
              Join[
                KeyValueMap["const DiracNet& " <> #2 <> "();"&, dMap],
                KeyValueMap["const NetVal& " <> #2 <> "();"&, lMap]],
              "\n"]],
        "DiracPool" -> dPoolCse, "LorPool" -> lPoolCse|>]];

(* ---- STAGE 5 log: the dedup join's summary (regen_check.sh greps the "[cse] sub-terms:" line) ---- *)
ntGenLogDedup[joined_Association, noDedup_] :=
  With[{subKeysLen = joined["subKeysLen"], netTraceRows = joined["netTraceRows"],
        nSub = joined["nSub"], nReused = joined["nReused"]},
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
      If[noDedup,
        " [NT_GEN_NO_DEDUP]",
        ""]]];

(* ---- STAGE 6: unit tables ---------------------------------------------------------------------
   The distinct traces are emitted as ONE flat table each, chunked by ntChunkDef and bin-packed into
   the -O0 unit TUs; the nets reference them by index. *)

(* DRESSED slot tables — INTERNED. Thousands of sub-terms share a handful of chains and draw their
   DSlotOpts from a tiny option pool, so per-sub-term literals would bloat the generator source ~20x.
   Emit the distinct chains (`chp`) and options (`optp`) ONCE, and each sub-term as INDICES (`sdchR`:
   its chain; `sdslR`: one option per slot); main() rebuilds sdch[k]/sdsl[k] in an O(nSub) loop.
   Non-dressed sub-terms carry an empty chain, keeping the sdch[k].empty() fast path. *)
ntGenDressedSlotTables[distinctSubs_, nSub_] :=
  Module[{chainStrs, combos, uChains, chainPos, uOpts, optPos, chainDefs, chainDecl, chainRefDefs,
          chainRefDecl, optDefs, optDecl, slotRefDefs, slotRefDecl},
    chainStrs = If[nSub === 0, {}, distinctSubs[[All, 3]]];
    (* per sub-term: its option-string LIST ({} for a non-slot sub-term) *)
    combos    = If[nSub === 0, {}, distinctSubs[[All, 4]]];
    uChains = DeleteDuplicates[chainStrs];
    chainPos = AssociationThread[uChains -> Range[Length[uChains]] - 1];
    uOpts = DeleteDuplicates[Flatten[combos]];
    optPos = AssociationThread[uOpts -> Range[Length[uOpts]] - 1];
    {chainDefs, chainDecl} = ntChunkDefs["chp", "std::vector<std::vector<DChainTok>>", {uChains}];
    {chainRefDefs, chainRefDecl} = ntChunkDefs["sdchR", "std::vector<int>", {ntIntStrs[chainPos /@ chainStrs]}];
    {optDefs, optDecl} = ntChunkDefs["optp", "std::vector<DSlotOpt>", {uOpts}];
    {slotRefDefs, slotRefDecl} = ntChunkDefs["sdslR", "std::vector<std::vector<int>>",
      {ntIntRow[optPos /@ #]& /@ combos}];
    ntStageResult["ntGenDressedSlotTables", {"Defs", "ChunkDecls"},
      <|"Defs" -> Join[chainDefs, chainRefDefs, optDefs, slotRefDefs],
        "ChunkDecls" -> chainDecl <> chainRefDecl <> optDecl <> slotRefDecl|>]];

(* the distinct-trace tables sdn (Dirac nets) and sln (Lorentz nets), plus the dressed slot tables *)
ntGenTraceTables[distinctSubs_, nSub_, hasDressed_] :=
  Module[{sdnDefs, sdnDecl, slnDefs, slnDecl, dressed},
    {sdnDefs, sdnDecl} = ntChunkDefs["sdn", "std::vector<DiracNet>", If[nSub === 0, {{}}, {distinctSubs[[All, 1]]}]];
    {slnDefs, slnDecl} = ntChunkDefs["sln", "std::vector<NetVal>", If[nSub === 0, {{}}, {distinctSubs[[All, 2]]}]];
    dressed = If[hasDressed, ntGenDressedSlotTables[distinctSubs, nSub], <|"Defs" -> {}, "ChunkDecls" -> ""|>];
    ntStageResult["ntGenTraceTables", {"Defs", "ChunkDecls"},
      <|"Defs" -> Join[sdnDefs, slnDefs, dressed["Defs"]],
        "ChunkDecls" -> sdnDecl <> slnDecl <> dressed["ChunkDecls"]|>]];

(* the shared decl header: net-builder + CSE-accessor forward declarations (the units that call the
   lc<k>()/dc<k>() accessors #include this — emitted ONCE here, not duplicated per unit). *)
ntGenDeclHeader[hasDressed_, cseDecls_String, chunkDecls_String] :=
  ntStageResult["ntGenDeclHeader", {"Decl"},
    <|"Decl" ->
      "// GENERATED by MakeNTKernel — do not edit. Numeric net-builder declarations.\n#pragma once\n" <> "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n#include <vector>\n" <>
        If[hasDressed,
          "#include \"numtracer/numeric/numeric_contract.hpp\"\n",
          ""
        ] <> "using namespace numtracer::network;\n" <>
        If[hasDressed,
          "using namespace numtracer::numeric;\n",
          ""
        ] <> cseDecls <> "\n" <>
(* the DISTINCT-trace tables (see ntGenDedupJoin): one flat builder each, which the nets index into. *)
        "std::vector<DiracNet> sdn0();\n" <> "std::vector<NetVal> sln0();\n" <>
        If[hasDressed,
          "std::vector<std::vector<DChainTok>> chp0();\n" <> "std::vector<int> sdchR0();\n" <>
          "std::vector<DSlotOpt> optp0();\n" <> "std::vector<std::vector<int>> sdslR0();\n",
          ""
        ] <>
(* chunk helpers of the oversized builders — the bin-packer may put a builder's helpers in a
   different unit than its assembler, so these must be declared here, not per-unit. *)
        chunkDecls|>];

(* Timed apart from the main() data tables: this half is dominated by ntChunkDef's dedup scans,
   that one by integer-to-text. *)
ntGenUnitSources[distinctSubs_, nSub_, hasDressed_, cseDefs_List, cseDecls_String, colChunkDefs_List, unitPre_String] :=
  Module[{traces, allDefs, nUnits, units, decl, ntT},
    ntT = First @ AbsoluteTiming[
      traces = ntGenTraceTables[distinctSubs, nSub, hasDressed];
      allDefs = Join[cseDefs, traces["Defs"], colChunkDefs];
(* Size-aware unit count (~$ntUnitChars per unit, 8..$ntUnitCap): a count-based rule would emit
   hundreds of tiny TUs (the CSE accessors are many but small), each re-parsing the shared header.
   Total[StringLength /@ ...] avoids materialising every def as one giant string just to measure it. *)
      nUnits = Min[Min[$ntUnitCap, Max[8, Ceiling[Total[StringLength /@ allDefs] / $ntUnitChars]]], Max[1, Length[allDefs]]];
      units =
        If[allDefs === {},
          {},
          (unitPre <> StringRiffle[#, "\n"] <> "\n")& /@ ntLptBinPack[allDefs, nUnits]];
      decl = ntGenDeclHeader[hasDressed, cseDecls, traces["ChunkDecls"]]["Decl"];];
    ntLog["[prof] unit tables + bin-pack: ", ntT, " s"];
    ntStageResult["ntGenUnitSources", {"Units", "Decl"}, <|"Units" -> units, "Decl" -> decl|>]];

(* ---- STAGE 7: main() --------------------------------------------------------------------------
   Each function returns main()'s text for its section; those owning big tables also return the
   table definitions ("Defs"), which the driver prepends to main(). *)

(* argv, the LorentzEnv, the component table and the distinct-trace tables *)
ntGenMainPrologue[ncomp_, nsInner_, hasDressed_] :=
  ntStageResult["ntGenMainPrologue", {"Text"},
    <|"Text" ->
      StringJoin[
        "int main(int argc, char** argv){\n",
        If[ntSingleQ[], "  numtracer::codegen::emit_precision() = numtracer::codegen::EmitPrecision::Single;\n", ""],
        "  std::string decor = \"static inline\"; std::string hns = \"" <> nsInner <> "\";\n",
        "  for(int a=1;a<argc;++a){ std::string s=argv[a]; if(s==\"-d\"&&a+1<argc) decor=argv[++a]; else if(s==\"-n\"&&a+1<argc) hns=argv[++a]; }\n",
        "  const int nsym = " <> ToString[ncomp["nsym"]] <> ";\n",
(* units is emitted BEFORE the LorentzEnv because the env binds both nsym and the unit groups;
   comp/atomDen and the trace entry points are then built through `env`. *)
        "  std::vector<std::vector<int>> units = {" <> StringRiffle[ntIntRow /@ ncomp["units"], ","] <> "};\n",
        "  LorentzEnv env(nsym, units);\n",
        "  std::vector<std::array<MPoly,4>> comp(" <> ToString[ncomp["maxBase"] + 1] <> ", {env.zero(),env.zero(),env.zero(),env.zero()});\n",
        (* component-table init: comp[base][mu] = <MPoly builder>, skipping structural zeros. *)
        KeyValueMap[
          Function[{base, comps},
            MapIndexed[
              Function[{s, mu},
                If[s === "env.zero()",
                  "",
                  "  comp[" <> ToString[base] <> "][" <> ToString[mu[[1]] - 1] <> "] = " <> s <> ";\n"]],
              comps]],
          ncomp["compCpp"]],
        "  std::vector<std::string> symNames = {" <> StringRiffle[("\"" <> # <> "\"")& /@ ncomp["symNamesCpp"], ","] <> "};\n",
(* the DISTINCT-trace tables (see ntGenDedupJoin). Flat, indexed by trace id: a plain trace has an
   empty chain (contract via sdn[k]), a dressed (structural) one uses sdch[k]/sdsl[k] via
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
          ""]]|>];

(* Per net: which traces it references (sidx), with what scalar (dsc) and, dressed flows only, which
   dressing monomial (sdr); sub-terms sharing a trace were already merged, and zero sums dropped, by
   ntGenDedupJoin. HASH-CONSED (ntHashConsRows): written out in full these tables dominate the main
   TU, and a multi-megabyte braced-init in the optimised main TU costs minutes. sidx is row-deduped;
   dsc and sdr on BOTH levels (distinct values, then distinct rows of value indices). The O(nets)
   runtime rebuild reproduces sidx/dsc/sdr exactly.
   Emitted names: <t>U = the distinct rows, <t>R = each net's row index into <t>U, <t>V = the
   distinct values a <t>U row indexes (ntSidxU, ntDscU, ... are the chunk functions building them). *)
ntGenMainTables[netTraceRows_, netDressRows_, netScalarRows_, hasDressed_] :=
  Module[{sidx = ntHashConsRows[netTraceRows], dsc = ntHashConsRows[netScalarRows, "Values"],
          sdr, sidxU, sidxR, dscU, dscR, sdrR},
    sidxU = ntGenBigTable["ntSidxU", "std::vector<std::vector<int>>", ntIntRow /@ sidx["Rows"], "std::vector<std::vector<int>> sidxU"];
    sidxR = ntGenBigTable["ntSidxR", "std::vector<int>", ntIntStrs[sidx["RowIdx"]], "std::vector<int> sidxR"];
    dscU = ntGenBigTable["ntDscU", "std::vector<std::vector<int>>", ntIntRow /@ dsc["Rows"], "std::vector<std::vector<int>> dscU"];
    dscR = ntGenBigTable["ntDscR", "std::vector<int>", ntIntStrs[dsc["RowIdx"]], "std::vector<int> dscR"];
    If[hasDressed,
      sdr = ntHashConsRows[netDressRows, "Values"];
      sdrR = ntGenBigTable["ntSdrR", "std::vector<int>", ntIntStrs[sdr["RowIdx"]], "std::vector<int> sdrR"]];
    ntStageResult["ntGenMainTables", {"Defs", "Text"},
      <|"Defs" -> If[hasDressed, {sidxU[[1]], sidxR[[1]], dscU[[1]], dscR[[1]], sdrR[[1]]}, {sidxU[[1]], sidxR[[1]], dscU[[1]], dscR[[1]]}],
        "Text" ->
          StringJoin[
            "  const size_t NNET = " <> ToString[Length[netTraceRows]] <> ";\n",
            sidxU[[2]],
            sidxR[[2]],
            "  std::vector<std::vector<int>> sidx(NNET);\n",
            "  for(size_t i=0;i<NNET;++i) sidx[i]=sidxU[sidxR[i]];\n",
            "  std::vector<Cx> dscV = {" <>
              StringRiffle[
                Function[s,
                    "Cx{" <> cppNum[Re[s]] <> "," <> cppNum[Im[s]] <> "}"
                  ] /@ dsc["Values"],
                ","
              ] <> "};\n",
            dscU[[2]],
            dscR[[2]],
            "  std::vector<std::vector<Cx>> dsc(NNET);\n",
            "  for(size_t i=0;i<NNET;++i){ const auto& r=dscU[dscR[i]]; dsc[i].reserve(r.size());\n",
            "    for(int k: r) dsc[i].push_back(dscV[k]); }\n",
(* the dressed phase-B fold routes each sub-term's scaled MPoly trace into its DPoly channel
   sdr[i][j] (empty monomial = undressed); sdrV and sdrU are small, so they stay inline. *)
            If[hasDressed,
              "  std::vector<DMono> sdrV = {" <>
                StringRiffle[ntIntRow /@ sdr["Values"], ","] <> "};\n" <>
              "  std::vector<std::vector<int>> sdrU = {" <>
                StringRiffle[ntIntRow /@ sdr["Rows"], ","] <> "};\n" <>
              sdrR[[2]] <>
              "  std::vector<std::vector<DMono>> sdr(NNET);\n" <>
              "  for(size_t i=0;i<NNET;++i){ const auto& r=sdrU[sdrR[i]]; sdr[i].reserve(r.size());\n" <>
              "    for(int k: r) sdr[i].push_back(sdrV[k]); }\n",
              ""]]|>]];

(* atom denominators, Matsubara-evenness probe, worker counts, the trace lambda and PHASE A *)
ntGenMainPhaseA[nSub_, nReused_, hasDressed_, mIdx_] :=
  ntStageResult["ntGenMainPhaseA", {"Text"},
    <|"Text" ->
      StringJoin[
(* the projector atom denominators are keyed by ATOM ID (e.inv/e.invS), not by position, and each
   id is filled idempotently — and the distinct traces cover every lnet that occurs. So scanning
   the deduped table gives the same atomDen as scanning every occurrence did. *)
        "  auto atomDen = env.collect_atom_denoms(sln, comp);\n",
        "  for(auto &a: atomDen) a = reduce_units(a, units);  // bare-loop k^2 -> monomial l1^2 -> cancels\n",
(* MATSUBARA EVENNESS, proven while contracting. If every trace and every atom denominator carries
   only EVEN powers of the Matsubara frequency, kernel(+w) == kernel(-w), and DiFfRG may collapse
   `kernel(+w) + kernel(-w)` to `2*kernel(w)`, halving the Matsubara-sum work.
   It must be proven HERE: DiFfRG's MakeKernel "MatsubaraEven" option tests the placeholder
   `body = 0.` (DiFfRG_compat.m), so it is trivially True for every flow and would silently drop the
   odd half of a non-even kernel. The test is SUFFICIENT only (cancelling odd terms read as odd),
   which is the safe direction: a false "odd" costs an optimisation, a false "even" is wrong physics.
   mIdx is the MPoly var index (0-based) of the Matsubara frequency, or -1 for "not finite-T". *)
        If[mIdx >= 0,
          "  // Matsubara evenness (see poly_even_in): every trace and every atom denominator must\n" <>
          "  // carry only even powers of var(" <> ToString[mIdx] <> "), the Matsubara frequency.\n" <>
          "  std::atomic<bool> ntMEven{true};\n" <>
          "  for(const auto &a: atomDen) if(!poly_even_in(a, " <> ToString[mIdx] <> ")) ntMEven.store(false, std::memory_order_relaxed);\n",
          ""],
        "  const bool ntprof = (std::getenv(\"NT_GEN_PROFILE\")!=nullptr);\n",
        "  unsigned workersA=std::thread::hardware_concurrency(); if(!workersA)workersA=4u;\n",
        "  if(const char* mw=std::getenv(\"NT_GEN_MAXW\")){int v=std::atoi(mw); if(v>0&&(unsigned)v<workersA)workersA=(unsigned)v;}\n",
(* SEPARATE worker count for phase B, whose memory profile differs from phase A's. Phase B runs `hw`
   concurrent RECOMPUTES of the uncached traces, which are the singletons (ntGenDedupJoin orders by
   descending refcount) and typically the heaviest ones; they land on top of the live window, so
   phase B can need FEWER workers than phase A. Defaults to phase A's count. *)
        "  unsigned workersB=workersA; if(const char* mb=std::getenv(\"NT_GEN_MAXW_B\")){int v=std::atoi(mb); if(v>0)workersB=(unsigned)v;}\n",
        "  const long NSUB = " <> ToString[nSub] <> ";\n",
(* how many traces are RESIDENT. Default: the reused ones (refCount >= 2), which the dedup ordering
   puts first; a singleton is contracted once whether cached or not, so caching it is pure RAM.
   DRESSED flows default to NSUB: their traces are individually small, and phase B is parallel over
   NETS, so leaving singletons to it lets one dominant net serialise thousands of contractions.
   NT_GEN_MEMO_MAX overrides either way (clamped to [0, NSUB]): lower it when memory is tight, raise
   it to put the singletons in phase A's balanced work list. *)
        "  long nCache = " <> ToString[If[hasDressed, nSub, nReused]] <> ";\n",
        "  if(const char* mm=std::getenv(\"NT_GEN_MEMO_MAX\")){ long v=std::atol(mm); if(v>=0) nCache=std::min<long>(v,NSUB); }\n",
(* Traces are dressing-stripped, so phase A caches plain MPoly on both paths; the DPoly is assembled
   only in the phase-B fold. *)
        "  auto trace=[&](int k)->MPoly{\n",
(* The parity probe sits on the trace lambda, not the trace table: with nCache == 0 the table is
   EMPTY and phase B recomputes through this lambda, so every distinct trace passes through here.
   DRESSED flows are excluded: their dressing monomials (sdr) are multiplied back in by the group fold
   and not seen here, and a dressing could carry an odd power. *)
        If[mIdx >= 0 && !hasDressed,
          "    MPoly tracePoly = env.numeric_value_netval(sdn[k], sln[k], comp, atomDen);\n" <>
          "    if(!poly_even_in(tracePoly, " <> ToString[mIdx] <> ")) ntMEven.store(false, std::memory_order_relaxed);\n" <>
          "    return tracePoly;\n",
          If[hasDressed,
            (* structural trace → plain MPoly: a non-slot sub-term contracts via sdn[k]/sln[k],
               a collected one via the dressing-free _mp variant (dressing stripped at codegen). *)
            "    return sdch[k].empty()\n" <> "      ? env.numeric_value_netval(sdn[k], sln[k], comp, atomDen)\n" <> "      : env.numeric_value_dressed_netval_mp(sdch[k], sdsl[k], sln[k], comp, atomDen);\n",
            "    return env.numeric_value_netval(sdn[k], sln[k], comp, atomDen);\n"]],
        "  };\n",
        (* PHASE A — contract each distinct trace once, parallel over a FLAT work list (numeric/trace_fold.hpp). *)
        "  auto tA=std::chrono::steady_clock::now();\n",
        "  std::vector<MPoly> traceTable = env.contract_traces<MPoly>(nCache, workersA, trace);\n",
        "  if(ntprof){ std::size_t tb=0; for(auto &p: traceTable) tb+=poly_bytes(p);\n",
        "    std::fprintf(stderr,\"[num] phase A: %ld distinct traces, %ld cached, table %.1f MB, %.1f s (W=%u)\\n\",\n",
        "      NSUB, nCache, tb/1048576.0, std::chrono::duration<double>(std::chrono::steady_clock::now()-tA).count(), workersA); }\n",
(* RELEASE PHASE A'S ARENA. Concurrent dense contractions transiently allocate GBs, which glibc keeps
   in its per-thread arenas, so phase B would start from phase A's RSS high-water mark rather than
   from the (much smaller) trace table. malloc_trim returns it to the OS in milliseconds. *)
        "#if defined(__GLIBC__)\n",
        "  { const double rssPre = ntRssMB(); malloc_trim(0);\n",
        "    if(ntprof) std::fprintf(stderr,\"[num] arena trim after phase A: RSS %.0f -> %.0f MB\\n\", rssPre, ntRssMB()); }\n",
        "#endif\n",
(* PHASE B is fused into the group/lowering loop (ntGenMainPhaseB), which needs `groups`, `colv` and
   `realOnly` declared first. `tB` starts here so the reported time covers it. *)
        "  auto tB=std::chrono::steady_clock::now();\n"]|>];

(* the groups table, the symbol environment, and PHASE B + group accumulation + lowering, FUSED into
   one streaming pass (numeric/trace_fold.hpp's fold_groups_streaming). `groups` partitions the nets,
   so a net's folded polynomial is needed by exactly one group and dies once that group has absorbed
   it; only a window of nets and group accumulators is ever live.
   ORDER INVARIANT: each group left-folds its members in group order, and the sink runs on the calling
   thread for gi = 0,1,2,... ascending; that fixes GlobalEnv intern order and the CSE instruction
   stream. The scale is spelled per branch so each keeps its exact expression (poly*constant vs
   DPoly scaleCx). With CrossTraceCSE (crossCSE) the sink feeds ONE shared CSE builder (FusedStream)
   instead of one independent program per group. *)
ntGenMainPhaseB[groups_, nNet_, realOnlyG_, hasDressed_, crossCSE_] :=
  With[{nGrp = Length[groups],
        grp = ntGenBigTable["ntGroups", "std::vector<std::vector<int>>", ntIntRow /@ groups, "std::vector<std::vector<int>> groups"]},
    ntStageResult["ntGenMainPhaseB", {"Defs", "Text"},
      <|"Defs" -> {grp[[1]]},
        "Text" ->
          StringJoin[
            grp[[2]],
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
                  If[TrueQ[#], "1", "0"]& /@ realOnlyG,
                  Table["0", {nGrp}]],
                ","
              ] <> "};\n",
            "  long netWindow = numtracer::numeric::net_window((long)sidx.size(), workersB);\n",
(* UNCONDITIONAL and FATAL (check_group_partition in numeric/trace_fold.hpp). `groups` must
   PARTITION the nets: a duplicate folds a net in twice, a gap silently drops one, and either is a
   wrong kernel with no other symptom. O(nNet), once per generation, so never gate it on a flag. *)
            "  numtracer::numeric::check_group_partition(groups, " <> ToString[nNet] <> ");\n",
(* The dressed fold reads the plain-MPoly trace table + the per-sub-term dressing monomials sdr and
   builds the per-net DPoly channel by channel (fold_groups_streaming_dressed). *)
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
            "  if(ntprof) std::fprintf(stderr,\"[num] phase B+lower: %d nets in %d groups, window %ld, %.1f s (W=%u)\\n\", " <> ToString[nNet] <> ", " <> ToString[nGrp] <> ", netWindow, std::chrono::duration<double>(std::chrono::steady_clock::now()-tB).count(), workersB);\n",
(* the trace table is dead once every group has folded; emission below only needs the lowered
   instruction streams, which are orders of magnitude smaller. *)
            "  { std::vector<MPoly> dead; traceTable.swap(dead); }\n"]|>]];

(* the fill formulas, the header preamble, the Matsubara verdict, and the trace bodies *)
ntGenMainEmission[varFill_, nsInner_, kernelNs_, fillArgSig_, complexQ_, hasDressed_, crossCSE_, mIdx_, nGrp_] :=
  ntStageResult["ntGenMainEmission", {"Text"},
    <|"Text" ->
      StringJoin[
(* Emission renders, dedups and writes every trace body serially (minutes on a multi-MB header);
   timed so the [num] trail covers the whole run. *)
        "  const auto tEmit = std::chrono::steady_clock::now();\n",
        "  FillFormulas fm;\n",
        "  fm.var = [](int id)->std::string{\n",
        Table["    if(id==" <> ToString[i - 1] <> ") return \"" <> varFill[[i]] <> "\";\n", {i, 1, Length[varFill]}],
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
   gcc _Complex builtins that device code silently computes as 0. A device consumer #defines
   NT_TRACE_COMPLEX (e.g. to cuda::std::complex<double>) BEFORE including the traces header; host
   consumers get std::complex<double>. The alias lives in the per-kernel namespace so multiple
   traces headers coexist in one TU. *)
        "  std::cout << \"#pragma once\\n#include <cmath>\\n" <>
          If[complexQ,
            "#include <complex>\\n",
            ""
          ] <> "namespace " <> kernelNs <> " { namespace \" << hns << \" {\\n\";\n",
        If[complexQ,
          "  std::cout << \"#ifndef NT_TRACE_COMPLEX\\n#define NT_TRACE_COMPLEX std::complex<" <> $ntRealT <> ">\\n#endif\\nusing nt_complex_t = NT_TRACE_COMPLEX;\\n\";\n",
          ""],
(* Unqualified fma/sqrt in a float header would otherwise bind the global double ::fma/::sqrt on
   the host, silently running that arithmetic in double. *)
        If[ntSingleQ[], "  std::cout << \"using std::fma;\\nusing std::sqrt;\\n\";\n", ""],
(* fill() emits negative powers, e.g. powr<-1>(1 - ...) in certain frames or situations:
   the helper must invert for N < 0, or every such factor silently becomes 1. *)
        "  std::cout << \"template<int N> \" << decor << \" " <> $ntRealT <> " powr(" <> $ntRealT <> " x){ " <> $ntRealT <> " r=" <> If[ntSingleQ[], "1.f", "1.0"] <> "; for(int i=0;i<(N<0?-N:N);++i) r*=x; return N<0?" <> If[ntSingleQ[], "1.f", "1.0"] <> "/r:r; }\\n\";\n",
        "  emit_env_layout(std::cout, genv);\n",
        "  std::cout << \"static inline constexpr int nenv = \" << genv.syms.size() << \";\\n\";\n",
(* The proven verdict, as a compile-time constant the kernel class picks up. Always emitted for a
   finite-T flow (rather than only when true) so the header records what was checked: a `false`
   here is a positive statement that the traces were tested and are not even, not a gap. *)
        If[mIdx >= 0,
          "  std::cout << \"// Matsubara evenness of the traces in var(" <> ToString[mIdx] <> "), proven from the monomial\\n\";\n" <>
          "  std::cout << \"// exponents at generation time (see poly_even_in). Consumed by the kernel class as\\n\";\n" <>
          "  std::cout << \"// DiFfRG's `matsubara_even` trait, which halves the Matsubara-sum evaluations.\\n\";\n" <>
          "  std::cout << \"static inline constexpr bool matsubara_even = \" << (" <>
            If[hasDressed, "false", "ntMEven.load(std::memory_order_relaxed)"] <>
            " ? \"true\" : \"false\") << \";\\n\";\n" <>
          "  if(ntprof) std::fprintf(stderr,\"[num] matsubara_even = %s\\n\", " <>
            If[hasDressed, "\"false (dressed flow: not checked)\"", "(ntMEven.load(std::memory_order_relaxed)?\"true\":\"false\")"] <> ");\n",
          ""],
        "  emit_fill(std::cout, genv, \"fill\", \"" <> fillArgSig <> "\", fm, decor);\n",
(* TRACE-BODY DEDUP: groups are keyed by dressing COEFFICIENT, finer than trace STRUCTURE, so many
   trN can have byte-identical bodies. Render each, key on the name-independent body, and emit a
   duplicate as a one-line forwarder `trN(f){ return trK(f); }` (K = first with that body); a pure
   source-size win. The forwarder's RETURN TYPE is read off the emitted signature (the token before
   " tr<i>("), since a complex return is spelled as the `nt_complex_t` alias. The decorator stays
   `decor`, so an out-of-lined canonical body does not drag `noinline` onto the forwarder.
   Fused (crossCSE): ONE trace_all(f, t[]) instead of nGrp trN(); the kernel reads tarr[i]. *)
        If[crossCSE,
          "  emit_cpp_fused(std::cout, fused, \"trace_all\", decor);\n",
          "  { std::unordered_map<std::string,std::string> seen; seen.reserve((size_t)" <> ToString[nGrp] <> ");\n" <>
          "    for(int i=0;i<" <> ToString[nGrp] <> ";++i){\n" <>
          "      const std::string nm = \"tr\"+std::to_string(i);\n" <>
          "      std::ostringstream os; emit_cpp(os, progs[i], nm, decor);\n" <>
          "      std::string s = os.str(); std::string body = s.substr(s.find('{'));\n" <>
          "      auto it = seen.find(body);\n" <>
          "      if(it==seen.end()){ seen.emplace(std::move(body), nm); std::cout << s; }\n" <>
          "      else { const std::string sig = s.substr(0, s.find(\" \"+nm+\"(\")); const std::string rt = sig.substr(sig.rfind(' ')+1);\n" <>
          "        std::cout << decor << \" \" << rt << \" \" << nm << \"(const " <> $ntRealT <> " *f) { return \" << it->second << \"(f); }\\n\"; } } }\n"
        ],
        "  std::cout << \"}} // namespace " <> kernelNs <> "::\" << hns << \"\\n\";\n",
        "  if(ntprof) std::fprintf(stderr,\"[num] emission: %.1f s\\n\", std::chrono::duration<double>(std::chrono::steady_clock::now()-tEmit).count());\n",
        "  return 0;\n}\n"]|>];

(* ---- the driver --------------------------------------------------------------------------------
   `mIdx` is the MPoly var index (0-based) of the Matsubara frequency, or -1 for "not a finite-T
   flow / unknown"; when >= 0 the generator proves Matsubara evenness while it contracts.
   NT_GEN_NO_DEDUP=1 turns the dedup join off (I4): every occurrence is its own trace, nothing is
   merged, dropped or cached. It is the escape hatch and the control for the equivalence test, which
   must compare kernel VALUES: dedup changes GlobalEnv interning order and so renumbers every sN. *)
emitNumericGenerator[coreNets_, restScalars_, colourNets_, groups_, ncomp_, nsInner_, fillArgSig_, kernelNs_:"numtracer_kernels", complexQ_:False, realOnlyG_ : {}, crossCSE_:False, mIdx_:-1] :=
  Module[{nNet = Length[coreNets], hasDressed = !FreeQ[coreNets, _ntDressedCore],
          noDedup = ntEnvFlag["NT_GEN_NO_DEDUP"], sub, col, preambles, cse, joined, unitSrc, main, ntT},
    sub = ntGenExpandSubTerms[coreNets, restScalars];
    col = ntGenColourTables[colourNets, nNet];
    preambles = ntGenPreambles[hasDressed, col["MainDecls"], col["ChunkDefs"] =!= {}];
    cse = ntGenNetCSE[sub["DiracNetPool"], sub["LorNetPool"], sub["DiracNetIds"], sub["LorNetIds"]];
    ntT = First @ AbsoluteTiming[
      joined =
        ntGenDedupJoin[
          sub["DiracNetIds"], sub["LorNetIds"], sub["SubScalars"], sub["ChainIds"], sub["SlotTupleIds"], sub["DressMonoIds"],
          cse["DiracPool"], cse["LorPool"], sub["ChainPool"], sub["SlotTuplePool"], sub["DressMonoPool"],
          hasDressed, noDedup];];
    ntLog["[prof] sub-term dedup join: ", ntT, " s"];
    sub = None;   (* the per-sub-term columns are dead after the join *)
    ntGenLogDedup[joined, noDedup];
    unitSrc = ntGenUnitSources[joined["distinctSubs"], joined["nSub"], hasDressed, cse["Defs"], cse["Decls"],
                col["ChunkDefs"], preambles["UnitPre"]];
(* main(). On a dense flow it is ~100% flat integer tables, so this timer measures integer-to-text
   throughput. The table definitions go in front of main(), in this order. *)
    ntT = First @ AbsoluteTiming[
      main =
        With[{tables = ntGenMainTables[joined["netTraceRows"], joined["netDressRows"], joined["netScalarRows"], hasDressed],
              phaseB = ntGenMainPhaseB[groups, nNet, realOnlyG, hasDressed, crossCSE]},
          StringJoin[
            tables["Defs"], col["ColRDef"], phaseB["Defs"],
            ntGenMainPrologue[ncomp, nsInner, hasDressed]["Text"],
            tables["Text"],
            ntGenMainPhaseA[joined["nSub"], joined["nReused"], hasDressed, mIdx]["Text"],
            col["MainText"],
            phaseB["Text"],
            ntGenMainEmission[ncomp["varFill"], nsInner, kernelNs, fillArgSig, complexQ, hasDressed, crossCSE, mIdx, Length[groups]]["Text"]]];];
    ntLog["[prof] main() data tables: ", ntT, " s"];
    {preambles["Pre"], unitSrc["Units"], unitSrc["Decl"], main, joined["nSub"]}];
