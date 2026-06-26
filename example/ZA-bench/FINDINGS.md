# Full-QCD AAqbq ZA — BEFORE rework (pre-existing generated artifacts)
NumTracer "Generate" (per-diagram trace + C++ combiner):
  kernel.hh   383 KB
  kernels.hh  217 KB   (273 trace fns, 6952 temporaries)
  generation  115.5 s  (from QCD_NumTracer/generation_times.m: NumTrace 0.2s + codegen 115.3s)
  compile     2.27 s   (g++ -O3 -march=native, standalone w/ shim), maxRSS 374 MB
  runtime     589.25 ns/eval  (3.26x FormTracer 180.73 ns)
FormTracer oracle:
  kernel.hh   19 KB    (1 collected kernel, 37 temporaries)
NT-vs-FORM agreement: ratio 0.82..1.31 pointwise (~20% integrated mismatch).
  => PRE-EXISTING: the QCD_NumTracer ZA derivation is unfinished/not yet faithful to the
     QCD_vacuum_147 FORM derivation. NOT introduced by the collection rework, which is an
     output-preserving regrouping. Benchmark gate = reworked-NT vs pre-rework-NT (~1e-12).

## NT-vs-FORM mismatch LOCALIZED (2026-06-14)
Compared the cos1-INTEGRATED result (not just pointwise). Sector isolation (zeroing vertex dressings):
  gauge only (gluon+ghost+tadpole) ... 3.8e-15  MATCH
  + AAqbq seagull (ZAAqbq1/2) ........ 2.4e-14  MATCH
  + quark loop ZAqbq1 only ........... 3.5e-15  MATCH
  + quark loop ZAqbq4 + ZAqbq7 ....... 7.7      WRONG  <-- the bug
Isolated ZAqbq4/ZAqbq7 contribution: NT/FORM ratio 0.85-0.92 (roughly constant, mildly
kinematics-dependent) => NumTracer UNDER-COUNTS the non-classical quark-gluon vertex structures
by ~12%; signature = a missing/mis-weighted tensor (cross-)term, not a global sign/factor.
Scope: PRE-EXISTING bug in QCD_NumTracer's handling of the AqbqDirect147 legs 4,7 in the ZA quark
loop (not covered by any in-tree test; only ZAqbq1 is). NOT my benchmark, NOT the dr-engine change
(ZA3 regen byte-identical). Independent of the size-blowup rework (the rework's nt-vs-nt_before
faithfulness gate holds regardless), but NT-vs-FORM correctness needs this fixed.

## Route-B GlobalCollect measured on the big QCD AAqbq ZA (2026-06-14)
                         kernel.hh   kernels.hh   total    compile   runtime         gc-vs-pertrace
  per-trace (before)     383 KB      217 KB       600 KB   2.27 s    593 ns (3.23x)  --
  GlobalCollect (Route B)  7 KB      266 KB       273 KB   1.36 s    776 ns (4.23x)  5.8e-6 (algebraically same)
=> kernel.hh 54x smaller, total 2.2x smaller, compile faster — BUT runtime REGRESSED 3.23x->4.23x
   and the single collected polynomial loses ~6 digits to cancellation.
WHY: flattening the dressings into ONE polynomial over {kinematic, dressing-atom} symbols INFLATES the
   monomial/op count (dressing x kinematic products expand) and wrecks conditioning. FORM does NOT
   flatten — it keeps a SUM of factored (dressing-combo x Horner-factored kinematic) terms (~500 ops).
   The 273-trace bloat moved into one slow 266 KB trace; op count went UP, not down.
LESSON: the FORM-compactness win needs STRUCTURED (factored) collection — keep per-diagram/per-dressing
   factoring and SHARE kinematic sub-polynomials + dressing atoms across diagrams via CSE (the
   CrossTraceCSE direction), NOT flatten. Route-B (flatten) is correct + shrinks kernel.hh but is a net
   runtime loss for this flow. Kept OPT-IN (GlobalCollect, default off) so nothing regresses.

## (a) SOLVED — structured GlobalCollect: fold colour numerically + group by dressing coefficient
The flatten (Route-B v1) and group-by-dressing-MONOMIAL approaches both LOST. The win is to group by
the diagram's FULL dressing coefficient with ALL colour folded numerically: the quark-loop colour
channels the legacy path kept as fundamental SINGLETONS now merge into their Feynman graph.
                       traces  total    compile  runtime        gen
  per-trace (before)   273     600 KB   2.27 s   593 ns (3.23x)  115 s
  GlobalCollect (this) 29      54 KB    0.95 s   205 ns (1.12x)  118 s
  FORM                 --      19 KB    --       183 ns (1.00x)  --
Output-preserving vs per-trace = 1.2e-11 (well-conditioned: per-group Horner kept). Generation NOT
blown up (compile 76s net-builders + reduce/lower 0.8s; same units + ulimit -v RAM cap as before).
Implementation: `"GlobalCollect"->True` just forces fundamental colour through the numeric SUNNet
(colv) instead of the et `_col` singleton path, so the EXISTING group-by-diagData emit merges them.
No flatten, no dr-symbols in the collector (the dr engine + inv_dress test remain as validated, unused
infra). Enable per flow by adding "GlobalCollect"->True to the emit[] MakeNTKernel call.
Residual 1.12x vs FORM = within-trace op count (FORM collects+Horner-factors a touch tighter); good enough.

## (b) qbq4/qbq7 bug — LOCALIZED to NumTrace of slashed-momentum vertex structures
AqbqDirect147 vertex tensor structures (dumped via TBGetBasisElement):
  struct 1 (ZAqbq1): i·gamma[rho]·transProj(-p2-p3)_{mu,rho}                    [MOMENTUM-FREE Dirac]  -> OK (1e-14)
  struct 2 (ZAqbq4): -gamma[mu']·gamma[rho]·transProj·(p2+p3)_{mu'}             [slashed quark mom]    -> WRONG
  struct 3 (ZAqbq7): (i/2)·[gamma·p2, gamma·p3]·gamma[rho]·transProj            [slashed quark mom]    -> WRONG
Evidence the bug is in the NumTrace FRONT-END (not reduction/collection/projector):
  - gc-inv == per-trace-inv to 1.2e-11, and BOTH miss FORM by 0.57 -> reduction/collection faithful.
  - struct 1 (momentum-free) matches FORM 1e-14; only the slashed-momentum structs 4,7 fail.
  - transProj appears in ALL three structs -> projector handling is fine (else struct 1 would fail too).
  => NumTrace mis-traces the slashed-QUARK-MOMENTUM (gamma·p2, gamma·p3) part of the quark-gluon vertex
     in the loop. This NumTrace path is UNTESTED in-tree (in-tree uses AqbqDirect1 = struct 1 only; Zq
     also uses only struct 1). Likely: loop-momentum routing into the vertex's vec[p2]/vec[p3], or the
     Dirac trace of the extra gamma·p factors. Next: symbolic diff of NumTrace vs FORM for one qbq4
     diagram, OR user confirmation of the expected trace of structs 4,7.

Discriminator: with ZAqbq4==ZAqbq7 (same function) the err stays 0.47 (not ~1e-14) => NOT a struct<->dressing
swap; structures 4 AND 7 are each genuinely mis-traced. Ruled out: reduction, collection, projector,
momentum resolution (resolveComponents is a general linear substitution), struct/dressing pairing.
=> The remaining cause is NumTrace's Dirac trace of the slashed-quark-momentum (gamma.p) vertex structures
in the loop — a path no in-tree test exercises. Fix needs a symbolic NumTrace-vs-FORM diff of one qbq4
diagram (or user confirmation of the expected trace of structs 4,7). NumTracer engine (Mathematica
NumTrace / Dirac), NOT the QCD_NumTracer setup.

## (b) ISOLATED to struct-7 (ZAqbq7) MASS terms — trace engine proven correct
Sector toggles (full flow, Mq!=0):  qbq4-only = 7.7e-15 (OK)   qbq7-only = 2.18 (WRONG)   Mq=0 = 1e-14 (OK).
=> the bug is ENTIRELY struct-7 (the (i/2)[gamma.p2, gamma.p3] gamma^rho transProj commutator vertex),
   and ONLY its mass-insertion terms.
Battery of identical FunKit tensor exprs traced by FormTracer`FormTrace AND NumTracer NumTrace+Dense/ET
(sample l1=1.3,cos1=0.4,p=0.9), ALL MATCH to 1e-5:
   tr[gmu gmu]=16; tr[(g.p)(g.q)]; tr[g.a g.b g.c g.d] (4-slashed); transProj^2=12;
   struct4 loop (massless gg)=101.95; struct4 M^2 term=-68.38; struct4 gM=0;
   struct1xstruct4 single-mass (both orderings)=34.94/-33.43; struct7xstruct7 M^2 loop=-8.36;
   struct1xstruct7 cross mass loop=9.10; full complex (i*s1 - s4) massive loop Re=4.28.
=> NumTracer's Dirac/Lorentz TRACE PRIMITIVE is correct for struct-7 in every isolated config. The i/2
   complex coeff is NOT the issue (qbq7-only is real: i*i=-1 overall; qbq4-only HAS imaginary cross terms
   yet is correct). Ruled out: elementary trace, projector, slashed momenta, complex handling, reduction,
   collection, struct-4.
=> The discrepancy is in the DERIVATION PIPELINE producing struct-7's mass terms: FromFunKit / PropParam /
   FRoute momentum routing of the quark leg momenta (p2,p3) into the ANTISYMMETRIC commutator, or the Mq
   dressing momentum-argument — a QCD_NumTracer setup-level diff vs the QCD.m FORM pipeline, NOT NumTrace.
   struct-4 is SYMMETRIC in p2<->p3 (vec[p2]+vec[p3]) so routing-insensitive; struct-7 is ANTISYMMETRIC
   ([gamma.p2,gamma.p3]) so a p2<->p3 routing/argument diff flips it. Next: symbolic diff of the FromFunKit
   struct-7 mass diagram vs FormTracer's FlowAA struct-7 term.
SIDE BUG found: the ET/Dense backends emit Mathematica `Complex(re,im)` (CForm) for a complex SCALAR
   coefficient (cppNum/CForm doesn't handle Complex) -> compile error. Only hit with a colourless complex
   probe; the production Generate path emits via CppForm and is unaffected. Worth fixing cppNum for robustness.

## (b) FIXED — root cause + one-line fix in Codegen.m
ROOT CAUSE (bisected, not the Dirac trace, not reduce/rebase/ibp): the inv-path NET CONSTRUCTION routing.
In mkGenerateKernel, the per-component routing chose the fast path `chunkLorInv[dt0]` vs `splitColourGroupsInv`
based on `dt0 = expandDiracComponent[comp["Factors"]]`. But expandDiracComponent only Selects TOP-LEVEL
gammas; struct-7's vertex is (i/2)[gamma.p2,gamma.p3] gamma^rho with the commutator gammas BURIED in a
nested Plus. expandDiracComponent built a broken (open) Dirac chain -> traced to dt0=0. dt0=0 is gamma-free,
so the routing took the fast path chunkLorInv[0] = {} -> fell back to the UNIT net konst(1) -> the entire
struct-7 Lorentz structure was DROPPED. (struct-4's gammas are top-level, only its vec[p2]+vec[p3] is a Plus,
so dt0 was a correct trace -> fast path fine. struct-1 has no momenta. Hence only struct-7 broke, and only
its mass terms made the dropped structure observable.)
EVIDENCE: net0 for the L7 struct7xstruct7 probe = konst(1) (1 PTerm, 0 elements); reduce() of it = the
constant 1; but splitColourGroupsInv[factors] (which Expand/distributes BEFORE tracing) builds the correct
rich 2-entry net. expandDiracComponent + brute-force Lorentz contraction matched FORM (4.5995), confirming
the Dirac trace itself is correct.
FIX (Codegen.m, the routing condition): also route to splitColourGroupsInv when a top-level Plus factor
contains a Dirac element: `|| !FreeQ[Cases[comp["Factors"], _Plus], _ntGamma|_ntGamma5|_ntDeltaDirac]`.
splitColourGroupsInv Expands first so the buried chain is traced correctly.
RESULT: full QCD AAqbq ZA now matches FormTracer to 1.3e-14 integrated (was 0.57). Runtime 1.11x FORM,
54 KB (GlobalCollect win preserved). L7 minimal repro: inv 0.25 -> 8.36 (correct). 29/29 ctest; the fix is a
no-op for flows without gammas-in-a-Plus (ZA3/ZA4/Zq/ZAqbq1 unaffected; ZA3 regen still matches FORM 1.17e-14).
SIDE BUG (separate, pre-existing): cppNum/CForm emits Mathematica `Complex(re,im)` for a complex SCALAR
coefficient (ET/Dense backends) -> invalid C++. Not hit by the production Generate path; worth hardening.
