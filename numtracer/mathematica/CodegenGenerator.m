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
(* The key columns arrive already interned (stage 1 of emitNumericGenerator), so nothing is hashed
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
   they share. Returns {pre, units, decl, main}. The program it emits is what actually contracts the
   traces and PRINTS the committed straight-line kernel header; nothing here contracts anything.

   The stages, in order; each consumes the one before it:

     1. SUB-TERM EXPANSION + INTERNING. Walk every net's cores, expand the dressed slot options into
        their Cartesian product, and intern the five key columns (Dirac net, Lorentz net, dressed
        chain, slot-option tuple, dressing-atom multiset) through ntMkIntern. Output: six ragged
        integer columns, one entry per (net, sub-term), plus the five pools they index.
     2. NET CSE. Net-builder strings that recur across sub-terms become shared `lc<k>()` / `dc<k>()`
        accessor functions, and the pools are rewritten to call them.
     3. DEDUP JOIN (ntGenDedupJoin, above). Merge sub-terms sharing a trace AND a dress channel,
        count references, and order the distinct traces by descending reference count. Output: the
        distinct-trace table and the per-net fold entries that index it.
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

   The `Tuples` in stage 1 runs LAST-SLOT-FASTEST, and stages 1, 3 and 5 all rely on that alignment. *)

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
   The ORDER is fixed by the left-to-right evaluation of the StringJoin arguments below and is
   load-bearing: it decides the order the table functions appear in the generator source. *)
    tableBag = Internal`Bag[];
(* Emit one big table: stuff its definition into the bag, and return the `<lhs> = <call>;` line that
   goes in main(). *)
    emitBigTable[nm_String, ret_String, rows_List, lhs_String] :=
      With[{t = ntBigTableFns[nm, ret, rows]},
        Internal`StuffBag[tableBag, t[[1]]];
        "  " <> lhs <> " = " <> t[[2]] <> ";\n"];
(* DRESSED nets (symbolic dressing collection): a core may be ntDressedCore[chainStr, slotsStr].
   Traces are dressing-stripped, so the trace table is PLAIN MPoly on both paths: a dressed sub-term's
   structural trace contracts via numeric_value_dressed_netval_mp, and its dressing rides the
   per-sub-term scalar (dsc) + monomial (sdr, atom ids), assembled into the per-net DPoly by the
   phase-B fold (fold_groups_streaming_dressed). Combinations differing only in dressing thus share
   ONE trace. The dressings themselves are kind-2 `dress` env leaves filled by fm.dress. *)
    hasDressed = !FreeQ[invNets, _ntDressedCore];
    (* shared wrapper templates so the net-builder strings compile. *)
    tmpl = "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal tproj(){ return projT(Mu,Nu,Lb,Inv); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal lproj(){ return projL(Mu,Nu,Lb,Inv); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int InvS> NetVal mproj(){ return projM(Mu,Nu,Lb,InvS); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int Inv,int InvS> NetVal eproj(){ return projE(Mu,Nu,Lb,Inv,InvS); }\n" <> "template<int Mu,int Nu> NetVal lmetric(){ return met(Mu,Nu); }\n" <> "template<int Lbl,int Base,int Mask> NetVal lvec(){ return vec(Lbl,Base); }\n" <> "template<int A,int B,int C,int D> NetVal leps(){ return epsilon(A,B,C,D); }\n" <> "inline NetVal konst(double c){ return NetVal{PTerm{Cx{c,0}, {}}}; }\n" <> "template<class L> struct litco;\n" <> "template<numtracer::Cx C> struct litco<numtracer::Lit<C>>{ static constexpr numtracer::Cx v=C; };\n" <> "template<class L> NetVal sc(NetVal x){ return scale(litco<L>::v, std::move(x)); }\n";
(* STAGE 1: sub-term expansion. Per net, a colour group is a SUM of sub-terms: invNets[i] = {core_b…}
   (a DiracNet literal for a gamma branch, a Lorentz NetVal for a gamma-free one), invRest[i] =
   {{rest_b, scal_b}…} parallel. Each branch yields {ds, ls, scal, dc, dl, dr} columns (usually one
   sub-term). A DRESSED branch expands the Cartesian product of its chain's slot options, each a
   triple {structStr, num, dr}: the structures (dl) form the trace key, ∏ num folds into the scalar,
   and the dress-atom multiset (dr) becomes the DPoly key. Non-dressed branches carry an empty dress key.

   The five key columns (Dirac net, Lorentz net, chain, structural option tuple, dress multiset) are
   interned HERE, once per distinct value, so the dedup join works on integers. Byte-identity:
   `ntMkIntern` assigns ids in first-appearance order along the walk nets -> cores -> combinations
   (`Tuples` last-slot-fastest), which is the order the join relies on (I1). *)
    {dsInt, dsGet} = ntMkIntern[];
    {lsInt, lsGet} = ntMkIntern[];
    {dcInt, dcGet} = ntMkIntern[];
    {dlInt, dlGet} = ntMkIntern[];
    {drInt, drGet} = ntMkIntern[];
(* The Cartesian expansion of one core's slot options depends on `slotOpts` ALONE, so it is memoised
   per distinct option set. Returns {dlIds, drIds, nums, n} from one call so the columns can never
   desync: `Tuples` and `Flatten[Outer[...]]` both vary the last slot fastest, and that alignment
   between a structure tuple and its numeric coefficient is load-bearing. *)
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
                          (* gamma branch: DiracNet + projector rest *)
                          StringStartsQ[nv, "DiracNet"],
                            {{dsInt[nv]}, {lsInt[lsStr]}, ntPackCx[{scal}], {dcInt["std::vector<DChainTok>{}"]}, {dlInt[{}]}, {drInt[{}]}},
                          (* gamma-free branch: whole net is the rest *)
                          True,
                            {{dsInt["DiracNet{}"]}, {lsInt[nv]}, ntPackCx[{scal}], {dcInt["std::vector<DChainTok>{}"]}, {dlInt[{}]}, {drInt[{}]}}]]],
                    {cores, rss[[All, 1]], rss[[All, 2]]}]]],
          {invNets, invRest}];]},
      ntLog["[prof] sub-term expansion: ", ntT, " s"]];
(* ---- colour-net table: chunk DEFINITIONS on the parallel -O0 units, assembler in the main TU ----
   The distinct colour nets are one `SUNNet{...}` literal each; on a large colour graph the table can
   dominate the main TU, which is the one serial -O1 compile. So the chunk functions ride the
   LPT-packed -O0 units like every other builder (colChunkDefs -> allDefs below), with external
   (non-static) linkage, and only forward decls + the assembler stay in the main TU. *)
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
(* colour-net chunk defs ride the units (see colChunkDefs above); they build SUNNet literals
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
(* STAGE 2: NET-LEVEL CSE. Dense projections emit the SAME net sub-term thousands of times (e.g. the
   σ^μν quark-gluon vertex: ~500x), ballooning the generator source and its -O0 compile. Hash-cons each
   recurring net term into a shared accessor `lc<k>()` / `dc<k>()` (a function-local `static const`,
   so it is also BUILT once at run time). Trivial/empty literals and unique terms stay inline. Runs on
   both paths; each use copies the shared NetVal/DiracNet, so the committed kernel is unchanged. *)
    cseDefs = {};
    cseDecls = "";
    Module[{dPool = dsGet[], lPool = lsGet[], lCnt, dCnt, lMap = <||>, dMap = <||>, li = 0, di = 0},
(* Counted by interned ID. `Range[0, n-1]` yields the counts in id order, i.e. first-appearance
   order, which fixes the lc<k>/dc<k> numbering. *)
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
(* The rewrite is a POOL substitution: every occurrence of a net term shares one pool entry, so this
   is O(#distinct). Injective, because an `lc<k>()` / `dc<k>()` accessor can never collide with a net
   literal, so the pool's first-appearance order — and every downstream id — is untouched. *)
      With[{
        ntT =
          First @
            AbsoluteTiming[
              dPoolCse = If[KeyExistsQ[dMap, #], dMap[#] <> "()", #]& /@ dPool;
              lPoolCse = If[KeyExistsQ[lMap, #], lMap[#] <> "()", #]& /@ lPool;]},
        ntLog["[prof] CSE ref-rewrite: ", ntT, " s"]];
      ntLog["[cse] net terms: lnet ", Total[lCnt], "->", Length[lMap], " distinct, dnet ", Total[dCnt], "->", Length[dMap], " distinct shared builders"]
    ];
(* STAGE 3: the dedup join (ntGenDedupJoin; design and invariants there).
   NT_GEN_NO_DEDUP=1 turns it off: every occurrence is its own trace, nothing is merged, dropped or
   cached. It is the escape hatch and the control for the equivalence test, which must compare kernel
   VALUES: dedup changes GlobalEnv interning order and so renumbers every sN, making a byte-diff
   meaningless. *)
    ntNoDedup = ntEnvFlag["NT_GEN_NO_DEDUP"];
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
(* Published for the later COMPILE step (CodegenKernel.m), which logs it next to the main-TU
   optimisation level. Set HERE, at the hand-off, so a stale value from a previous flow never leaks
   into the next one. *)
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
(* STAGE 4: UNIT TABLES. The distinct traces are emitted as ONE flat table each, chunked by
   ntChunkDef and bin-packed into the -O0 unit TUs; the nets reference them by index. Timed apart
   from the main() data tables: this half is dominated by ntChunkDef's dedup scans, that one by
   integer-to-text. *)
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
(* DRESSED slot tables — INTERNED. Thousands of sub-terms share a handful of chains and draw their
   DSlotOpts from a tiny option pool, so per-sub-term literals would bloat the generator source ~20x.
   Emit the distinct chains (`chp`) and options (`optp`) ONCE, and each sub-term as INDICES (`sdchR`:
   its chain; `sdslR`: one option per slot); main() rebuilds sdch[k]/sdsl[k] in an O(nSub) loop.
   Non-dressed sub-terms carry an empty chain, keeping the sdch[k].empty() fast path. *)
    {chpDefs, chpCDecl, sdchrDefs, sdchrCDecl, optpDefs, optpCDecl, sdslrDefs, sdslrCDecl} =
      If[hasDressed,
        Module[{chainStrs, combos, uChains, chainPos, uOpts, optPos, chainDefs, chainDecl, chainRefDefs, chainRefDecl, optDefs, optDecl, slotRefDefs, slotRefDecl},
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
          {chainDefs, chainDecl, chainRefDefs, chainRefDecl, optDefs, optDecl, slotRefDefs, slotRefDecl}],
        {{}, "", {}, "", {}, "", {}, ""}];
    chunkDecls = sdnCDecl <> slnCDecl <> chpCDecl <> sdchrCDecl <> optpCDecl <> sdslrCDecl;
    allDefs = Join[cseDefs, sdnDefs, slnDefs, chpDefs, sdchrDefs, optpDefs, sdslrDefs, colChunkDefs];
(* Size-aware unit count (~$ntUnitChars per unit, 8..$ntUnitCap): a count-based rule would emit
   hundreds of tiny TUs (the CSE accessors are many but small), each re-parsing the shared header.
   Defs are packed by LPT (largest first, into the currently-smallest unit), since def sizes are
   skewed; ntChunkDef bounds any single def to ~$ntDefChunk, so the target is reachable. The cap must
   stay above what the size target asks for, or units quietly grow back past it.
   Total[StringLength /@ ...] avoids materialising every def as one giant string just to measure it. *)
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
(* the DISTINCT-trace tables (see ntGenDedupJoin): one flat builder each, which the nets index into. *)
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
(* STAGES 5-6: the main TU. On a dense flow it is ~100% flat integer tables, so this timer measures
   integer-to-text throughput. *)
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
(* the DISTINCT-trace tables (see ntGenDedupJoin). Flat, indexed by trace id: a plain
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
(* per net: which traces it references (sidx), and with what scalar (dsc); sub-terms sharing a trace
   were already merged, and zero sums dropped, by ntGenDedupJoin.
   HASH-CONSED: written out in full these tables dominate the main TU, and a multi-megabyte
   braced-init in the optimised main TU costs minutes. Both are extremely redundant, so dedupe rows
   for sidx, and for dsc BOTH levels (distinct values, then distinct rows of value-indices; rows alone
   differ while their entries repeat). The O(nets) runtime rebuild reproduces sidx/dsc exactly. *)
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
(* dressed flows only: the per-sub-term dressing monomials sdr[i][j], deduped on BOTH levels like dsc
   (few distinct monomials sdrV, few distinct index rows sdrU). The dressed phase-B fold routes each
   sub-term's scaled MPoly trace into its DPoly channel sdr[i][j] (empty monomial = undressed). *)
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
   only EVEN powers of the Matsubara frequency, kernel(+w) == kernel(-w), and DiFfRG may collapse
   `kernel(+w) + kernel(-w)` to `2*kernel(w)`, halving the Matsubara-sum work.
   It must be proven HERE: DiFfRG's MakeKernel "MatsubaraEven" option tests the placeholder
   `body = 0.` (DiFfRG_compat.m), so it is trivially True for every flow and would silently drop the
   odd half of a non-even kernel. The test is SUFFICIENT only (cancelling odd terms read as odd),
   which is the safe direction: a false "odd" costs an optimisation, a false "even" is wrong physics. *)
            If[mIdx >= 0,
              "  // Matsubara evenness (see poly_even_in): every trace and every atom denominator must\n" <>
              "  // carry only even powers of var(" <> str[mIdx] <> "), the Matsubara frequency.\n" <>
              "  std::atomic<bool> ntMEven{true};\n" <>
              "  for(const auto &a: atomDen) if(!poly_even_in(a, " <> str[mIdx] <> ")) ntMEven.store(false, std::memory_order_relaxed);\n",
              ""],
            "  const bool ntprof = (std::getenv(\"NT_GEN_PROFILE\")!=nullptr);\n",
            "  unsigned workersA=std::thread::hardware_concurrency(); if(!workersA)workersA=4u;\n",
            "  if(const char* mw=std::getenv(\"NT_GEN_MAXW\")){int v=std::atoi(mw); if(v>0&&(unsigned)v<workersA)workersA=(unsigned)v;}\n",
(* SEPARATE worker count for phase B, whose memory profile differs from phase A's. Phase B runs `hw`
   concurrent RECOMPUTES of the uncached traces, which are the singletons (ntGenDedupJoin orders by
   descending refcount) and typically the heaviest ones; they land on top of the live window, so
   phase B can need FEWER workers than phase A. Defaults to phase A's count. *)
            "  unsigned workersB=workersA; if(const char* mb=std::getenv(\"NT_GEN_MAXW_B\")){int v=std::atoi(mb); if(v>0)workersB=(unsigned)v;}\n",
            With[{
              (* traces are dressing-stripped, so phase A caches plain MPoly on both paths; the
                 DPoly is assembled only in the phase-B fold. *)
              PT = "MPoly"},
              StringJoin[
                {
                  "  const long NSUB = " <> str[nSub] <> ";\n",
(* how many traces are RESIDENT. Default: the reused ones (refCount >= 2), which the dedup ordering
   puts first; a singleton is contracted once whether cached or not, so caching it is pure RAM.
   DRESSED flows default to NSUB: their traces are individually small, and phase B is parallel over
   NETS, so leaving singletons to it lets one dominant net serialise thousands of contractions.
   NT_GEN_MEMO_MAX overrides either way (clamped to [0, NSUB]): lower it when memory is tight, raise
   it to put the singletons in phase A's balanced work list. *)
                  "  long nCache = " <> str[If[hasDressed, nSub, nReused]] <> ";\n",
                  "  if(const char* mm=std::getenv(\"NT_GEN_MEMO_MAX\")){ long v=std::atol(mm); if(v>=0) nCache=std::min<long>(v,NSUB); }\n",
                  "  auto trace=[&](int k)->" <> PT <> "{\n",
(* The parity probe sits on the trace lambda, not the trace table: with nCache == 0 the table is
   EMPTY and phase B recomputes through this lambda, so every distinct trace passes through here.
   DRESSED flows are excluded: their dressing monomials (sdr) are multiplied back in by the group fold
   and not seen here, and a dressing could carry an odd power. *)
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
(* RELEASE PHASE A'S ARENA. Concurrent dense contractions transiently allocate GBs, which glibc keeps
   in its per-thread arenas, so phase B would start from phase A's RSS high-water mark rather than
   from the (much smaller) trace table. malloc_trim returns it to the OS in milliseconds. *)
                  "#if defined(__GLIBC__)\n",
                  "  { const double rssPre = ntRssMB(); malloc_trim(0);\n",
                  "    if(ntprof) std::fprintf(stderr,\"[num] arena trim after phase A: RSS %.0f -> %.0f MB\\n\", rssPre, ntRssMB()); }\n",
                  "#endif\n",
(* PHASE B is fused into the group/lowering loop below (fold_groups_streaming), which needs
   `groups`, `colv` and `realOnly` declared first. `tB` starts here so the reported time covers it. *)
                  "  auto tB=std::chrono::steady_clock::now();\n"}]],
(* colnets HASH-CONSED: nets differing only in their Lorentz/Dirac part share a colour net, so
   sun_value_cx runs once per DISTINCT net and colR maps each net to it. The distinct table is built
   by ntColNets() in the preamble (see colMainDecls), off the main() body. *)
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
(* PHASE B + group accumulation + lowering, FUSED into one streaming pass (numeric/trace_fold.hpp's
   fold_groups_streaming). `groups` partitions the nets, so a net's folded polynomial is needed by
   exactly one group and dies once that group has absorbed it; only a window of nets and group
   accumulators is ever live.
   ORDER INVARIANT: each group left-folds its members in group order, and the sink runs on the calling
   thread for gi = 0,1,2,... ascending; that fixes GlobalEnv intern order and the CSE instruction
   stream. The scale is spelled per branch so each keeps its exact expression (poly*constant vs
   DPoly scaleCx).
   With CrossTraceCSE (crossCSE) the sink feeds ONE shared CSE builder (FusedStream) instead of one
   independent program per group. *)
            "  long netWindow = numtracer::numeric::net_window((long)sidx.size(), workersB);\n",
(* UNCONDITIONAL and FATAL (check_group_partition in numeric/trace_fold.hpp). `groups` must
   PARTITION the nets: a duplicate folds a net in twice, a gap silently drops one, and either is a
   wrong kernel with no other symptom. O(nNet), once per generation, so never gate it on a flag. *)
            "  numtracer::numeric::check_group_partition(groups, " <> str[nNet] <> ");\n",
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
            "  if(ntprof) std::fprintf(stderr,\"[num] phase B+lower: %d nets in %d groups, window %ld, %.1f s (W=%u)\\n\", " <> str[nNet] <> ", " <> str[nGrp] <> ", netWindow, std::chrono::duration<double>(std::chrono::steady_clock::now()-tB).count(), workersB);\n",
(* the trace table is dead once every group has folded; emission below only needs the lowered
   instruction streams, which are orders of magnitude smaller. *)
            "  { std::vector<MPoly> dead; traceTable.swap(dead); }\n",
(* Emission renders, dedups and writes every trace body serially (minutes on a multi-MB header);
   timed so the [num] trail covers the whole run. *)
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
   gcc _Complex builtins that device code silently computes as 0. A device consumer #defines
   NT_TRACE_COMPLEX (e.g. to cuda::std::complex<double>) BEFORE including the traces header; host
   consumers get std::complex<double>. The alias lives in the per-kernel namespace so multiple
   traces headers coexist in one TU. *)
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
(* TRACE-BODY DEDUP: groups are keyed by dressing COEFFICIENT, finer than trace STRUCTURE, so many
   trN can have byte-identical bodies. Render each, key on the name-independent body, and emit a
   duplicate as a one-line forwarder `trN(f){ return trK(f); }` (K = first with that body); a pure
   source-size win. The forwarder's RETURN TYPE is read off the emitted signature (the token before
   " tr<i>("), since a complex return is spelled as the `nt_complex_t` alias. The decorator stays
   `decor`, so an out-of-lined canonical body does not drag `noinline` onto the forwarder.
   Fused (crossCSE): ONE trace_all(f, t[]) instead of nGrp trN(); the kernel reads tarr[i]. *)
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
