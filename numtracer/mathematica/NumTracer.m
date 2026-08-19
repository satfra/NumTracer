(* ::Package:: *)

(* NumTracer` — Mathematica front-end that turns a traced tensor network into a
   C++ kernel driving the compile-time `et` engine (include/numtracer/).

   Two public entry points:
     NumTrace[net, opts]     builds an intermediate-expression tree (an NTKernel)
                             from a network written in the DSL heads below;
     MakeNTKernel[ntk, ...]  serialises that tree to C++ via the numeric
                             matrix-product backend (the sole generation path).

   The network is `Σ_terms coeff(dressings…) × trace(momenta)`: the trace is a
   contraction of the tensor heads (contracted numerically to an MPoly, then
   Horner-lowered to a real straight-line kernel), the coeff is flat C++. This is exactly the split
   tests/refshim/flow_ym_zc.hpp demonstrates by hand — we generate it. *)

BeginPackage["NumTracer`"];

Unprotect["NumTracer`*"];

ClearAll["NumTracer`*"];

(* ---- public entry points ---- *)

NumTrace::usage = "NumTrace[net, \"Frame\"->frame, \"Args\"->{...}] analyses a DSL tensor network into an NTKernel[<|...|>]: a list of terms, each a scalar coefficient times one or more independent contraction components, plus the env-id layout and frame needed to emit code.";

MakeNTKernel::usage = "MakeNTKernel[ntk, genFile, kernelFile, tracesFile, \"Name\"->\"X\", \"Namespace\"->ns, \"Dressings\"->{...}] serialises an NTKernel to C++ via the numeric matrix-product backend: a build-time generator program at `genFile`, run to produce the committed straight-line traces header `tracesFile`, and the kernel header `kernelFile` that fills the fundamental symbols and calls them.";

NTKernel::usage = "NTKernel[assoc] is the analysed intermediate-expression tree produced by NumTrace and consumed by MakeNTKernel.";

(* Attached to the package's primary entry point rather than to a symbol of its own: the file that
   failed may be the one that would have defined any more specific symbol. See ntLoadPart below. *)
NumTrace::loadsyntax = "NumTracer FAILED TO LOAD: `1` has a syntax error (see the Syntax:: message just above for the line). Mathematica does not treat this as a failure — Get returns Null, every definition before the malformed expression is installed and every one after it silently is not — so loading would otherwise continue and the package would be SILENTLY HALF-DEFINED. The usual symptom is far away and unrecognisable: a missing Options[] surfacing as SetOptions::optnf, a generator that emits nothing yet exits 0. An unbalanced bracket swallows everything to the end of the file, so the reported line is where parsing gave up, NOT where the error is — look for the unclosed bracket earlier. Aborting instead.";

FromFunKit::usage = "FromFunKit[flow, \"FlavourGroup\"->n] rewrites a FunKit traced flow (after // dressingRules) into the NumTracer DSL; scalar products become ntSP, resolved by the frame. SU(N) group tokens get their rank baked into the ntSUN* heads — colour from Global`Nc, the isospin group from \"FlavourGroup\" (default Global`Nf or 2).";

MakeNTKernelDiFfRG::usage = "MakeNTKernelDiFfRG[ntk, \"Name\"->\"ZA\", \"Integrator\"->\"Integrator_p2_1ang\", \"Parameters\"->kernelParameterList, \"IntegrationVariables\"->{\"l1\",\"cos1\"}, \"AngleDefs\"->{...}] scaffolds the DiFfRG flows/<Name>/ plumbing (MakeKernel on a placeholder body) and overwrites kernel.hh/kernels.hh with the NumTracer-traced numeric kernel. Dressings and their interpolator type are auto-derived from Parameters; the Regulator alias and the DiFfRG emission constants (namespace, includes, support API) are baked in. A loop-independent term flat-added to the integral (DiFfRG's constExpr) can be passed either positionally — MakeNTKernelDiFfRG[ntk, constExpr, opts] — or as \"Constant\"->expr; it populates constant(p,k,dressings) and is a plain expression (e.g. ZA[p]), not an NTKernel. Pair with UpdateNTFlows. Requires DiFfRG to be loaded.";

UpdateNTFlows::usage = "UpdateNTFlows[name] runs DiFfRG's UpdateFlows[name] then idempotently patches flows/CMakeLists.txt (find_package(NumTracer) + UNITY_BUILD OFF). Bundling the two makes the CMake patch atomic, so it can never be left un-applied after an UpdateFlows regenerates the file. Requires DiFfRG to be loaded.";

SetNumTracerThreads::usage = "SetNumTracerThreads[n] caps the code-generator to n worker threads for BOTH generation phases; SetNumTracerThreads[nA, nB] sets them separately. Phase A is the parallel Dirac-trace contraction (peak RSS ~nA x the per-contraction working set). Phase B is the per-net fold AND the CSE/Horner LOWERING of the kernel. Two consequences for nB: (1) each phase-B worker holds one full (possibly large) contraction, so on the dense quark flows a large nB drives peak memory; but (2) the lowering is compute-heavy and parallel, so a SMALL nB badly slows compute-bound flows with big kernels (e.g. the pure-gauge ZA4 / ZA4-full-basis: its phase-B+lower is ~600 s at nB=2 vs tens of seconds unthrottled). Therefore throttle ONLY around the RAM-heavy flows (the 4-point two-quark-two-gluon vertices) and keep nB high (or unset) everywhere else — do NOT set a low global cap at the top of a multi-flow script. The counts are exported as the NT_GEN_MAXW / NT_GEN_MAXW_B environment variables, read by the generator at run time (inherited through Run[]); call it any time BEFORE code generation. Under \"Offline\" -> True there is no Run[] to inherit them — the generator is launched later by the `numtrace` CMake target — so the values in force at emit time are additionally RECORDED into the flow's numtrace.json (\"maxw\"/\"maxw_b\", 0 = unset) and re-applied per flow by the build, where they lower (never raise) the build's own -jN. Returns {nA, nB}. SetNumTracerThreads[] reports the current setting. NT_GEN_MAXW only ever LOWERS the count below hardware concurrency; it cannot raise it.";

GetNumTracerThreads::usage = "GetNumTracerThreads[] returns {nA, nB} = the currently configured {phase-A, phase-B} generator thread caps (Automatic when unset, i.e. hardware concurrency). See SetNumTracerThreads.";

(* ---- DSL tensor heads (inert tags; first arg of the momentum-bearing ones is
        a momentum/4-vector symbol, the rest are contraction index labels) ---- *)

ntMetric::usage = "ntMetric[mu, nu] — Euclidean Lorentz metric delta_{mu nu}.";

ntVec::usage = "ntVec[q, mu] — tensor leg q_mu of momentum q (mu a symbolic Lorentz label). ntVec[q, i] with a literal INTEGER i (0-based, 0 = temporal/Matsubara) is the scalar component q_i, resolved by the frame like ntSP — used for finite-T objects such as ntVec[p, 0] = pi T.";

ntEpsFund::usage = "ntEpsFund[N, i1, ..., iN] -- the SU(N) FUNDAMENTAL Levi-Civita (totally antisymmetric invariant), carrying exactly N indices (colour SU(3): 3; isospin SU(2): 2). Comes from FunKit's epsFundCol/epsFundFlav. Never reaches the C++ engine: NumTrace contracts epsilon PAIRS into Kronecker deltas (eps.eps = k! det(delta)) via expandFundEps, since a lone epsilon is not an SU(N) invariant. The ADJOINT epsilon is separate and SU(2)-only -- there eps^abc = f^abc exactly, so FromFunKit rewrites it to ntSUNf[2,...]; that identification does NOT generalise to SU(N).";

ntUnitVec::usage = "ntUnitVec[i] — the constant unit basis 4-vector e_i (i = 0..3, 0 = temporal/Matsubara). Not written by hand: NumTrace introduces it when it rewrites a FIXED-component Lorentz index (gamma^0 and friends, the finite-T 3+1 split used by the four-quark Fierz bases) into a contraction with e_i, and injects its components into the frame. A fixed-component gamma is therefore emitted as an ordinary slash.";

ntSpatialVec::usage = "ntSpatialVec[q] — the SPATIAL part of momentum q as a momentum in its own right: components {0, q_1, q_2, q_3} (slot 0 = temporal/Matsubara). Not written by hand: FromFunKit introduces it for FormTracer's finite-T acronym vecs[q, mu], which becomes ntVec[ntSpatialVec[q], mu]. NumTrace pushes it through sums (it is linear) and injects the resulting components into the frame, so downstream it is an ORDINARY momentum leaf — in particular a spatial slash vecs[q,mu] gamma[mu,d1,d2] is emitted as an ordinary dslash and frameMask prunes its zero temporal component. Contrast ntSPS, the spatial scalar product, which is a scalar coefficient and needs no leaf.";

ntTransProj::usage = "ntTransProj[q, mu, nu] — transverse projector P_{mu nu}(q) = delta - q_mu q_nu/q^2 (valid at finite T).";

ntLongProj::usage = "ntLongProj[q, mu, nu] — longitudinal projector q_mu q_nu/q^2.";

ntElectricProj::usage = "ntElectricProj[q, mu, nu] — finite-T electric (time-like-transverse) projector P_E = P_T - P_M. Component 0 is the heat-bath/temporal direction.";

ntMagneticProj::usage = "ntMagneticProj[q, mu, nu] — finite-T magnetic (spatial-transverse) projector P_M_{ij} = delta_{ij} - q_i q_j/|q_vec|^2 (i,j spatial; temporal rows/cols vanish).";

(* ---- SU(N) group heads (colour, flavour/isospin, or any further group): one
        N-parameterized family. The leading argument N is the SU(N) rank, so several
        groups of different (or even equal) rank coexist in one network — the engine keeps
        them apart by their disjoint contraction ids and contracts each rank separately. ---- *)

ntSUNf::usage = "ntSUNf[N, a, b, c] — SU(N) structure constant f^{abc} (adjoint indices). The leading argument N is the SU(N) rank of the group the indices belong to.";

ntSUNDeltaAdj::usage = "ntSUNDeltaAdj[N, a, b] — SU(N) adjoint Kronecker delta^{ab} (rank N).";

ntSUNT::usage = "ntSUNT[N, a, i, j] — SU(N) fundamental generator (T^a)_{ij} (adjoint a, fundamental i,j; rank N).";

ntSUNDeltaFund::usage = "ntSUNDeltaFund[N, i, j] — SU(N) fundamental Kronecker delta_{ij} (rank N).";

ntSUNDiagFund::usage = "ntSUNDiagFund[N, i, j, spec] — a fundamental Kronecker delta_{ij} (rank N) carrying a PER-COMPONENT dressing. `spec` is a rules list {c1 -> expr1, c2 -> expr2, ..., Default -> defExpr}: `ci` are 1-based component indices (1..N) and each `expri` is a COMPLETE scalar dressing expression, kinematics included — e.g. {1 -> Zu[scale], 2 -> Zd[scale]}. A Default -> defExpr rule dresses every unnamed component, and components with neither a rule nor a Default are DROPPED (contribute nothing). Folds (via sun_value_dressed) to Σ_i c_i expr_i instead of the flavour-blind δ — this dresses e.g. the u- and d-quark differently WITHIN the SU(N) trace, no per-flavour diagram split. Each component may use its own scale. Kinematic heads (ntVec/ntSP/ntSPS) inside an expr are frame-resolved at kernel generation; a component that is a pure PROJECTOR rather than a carrier of physics is simply `c -> 1`.";

ntSUNDiagAdj::usage = "ntSUNDiagAdj[N, a, b, spec] — an adjoint Kronecker delta^{ab} (rank N) carrying a PER-COMPONENT adjoint dressing. `spec` is a rules list {c1 -> expr1, ..., Default -> defExpr} over the N^2-1 adjoint components, each `expri` a COMPLETE scalar dressing expression. This is how a gluon condensed along the Cartan is written (ntSUNDiagAdj[3, a, b, {3 -> A03[scale], 8 -> A08[scale]}]), the other 6 colours dropping out with no dead terms. Pinning a FIXED adjoint index uses the same head with `c -> 1`, since SUNFac has no pinned-index kind of its own.";


(* ANONYMOUS blanks, not named patterns. These two definitions live in the PUBLIC context (above
   Begin["`Private`"]), so every pattern NAME they use is created as an exported symbol — `n`, `i`,
   `j`, `a`, `b`, `spec`, `scale` — which `Protect["NumTracer`*"]` at the foot of this file then
   makes read-only. Any consumer that legitimately assigns to `Global`b` or `Global`n` after loading
   NumTracer then hits `Set::wrsym: Symbol b is Protected` and is SILENTLY refused, because the
   exported symbol shadows theirs on $ContextPath. That is not hypothetical: it is what stopped
   tests/gen/gen_fierz_ortho_numeric.wls from loading its basis fixture, and because that gate
   correctly refuses to pass vacuously it failed the run — and took regen_check.sh's whole test phase
   with it. The bodies below reference none of these names, so blanks are a pure deletion. *)
(* The 5-argument (N, i, j, spec, scale) spelling is gone. Reject it AT CONSTRUCTION rather than
   letting it reach the emitter: colFacG only runs at kernel-generation time, so a stale call would
   otherwise survive NumTrace and fail much later, far from the line that wrote it. *)
ntSUNDiagFund[_, _, _, _, _] := (
    Print["[NumTracer] ERROR: ntSUNDiagFund no longer takes a separate `scale`. Apply each dressing ",
      "to its own kinematics in the spec:\n",
      "    old:  ntSUNDiagFund[N, i, j, {1 -> Zu, 2 -> Zd}, scale]\n",
      "    new:  ntSUNDiagFund[N, i, j, {1 -> Zu[scale], 2 -> Zd[scale]}]\n",
      "  (a pure projector component is now `c -> 1`, not `c -> ntUnitDressing`.)"];
    Abort[]);

ntSUNDiagAdj[_, _, _, _, _] := (
    Print["[NumTracer] ERROR: ntSUNDiagAdj no longer takes a separate `scale`. Apply each dressing ",
      "to its own kinematics in the spec:\n",
      "    old:  ntSUNDiagAdj[N, a, b, {3 -> A03, 8 -> A08}, scale]\n",
      "    new:  ntSUNDiagAdj[N, a, b, {3 -> A03[scale], 8 -> A08[scale]}]\n",
      "  (a pure projector component is now `c -> 1`, not `c -> ntUnitDressing`.)"];
    Abort[]);

ntSP::usage = "ntSP[q1, q2] — Lorentz scalar product q1.q2 (a scalar coefficient).";

ntSPS::usage = "ntSPS[q1, q2] — spatial (finite-T) scalar product q1_vec.q2_vec = ntSP[q1,q2] - q1_0 q2_0 (a scalar coefficient).";

ntDress::usage = "ntDress[head, args...] — an opaque scalar dressing emitted verbatim as head(args).";

(* ---- Dirac (spinor) sector heads. din/dout are spinor axis labels; mu a Lorentz label.
        A slashed momentum is ntGamma[mu, din, dout] contracted with ntVec[q, mu]. ---- *)

ntGamma::usage = "ntGamma[mu, din, dout] — Dirac gamma^mu carrying a Lorentz axis mu and spinor axes din,dout (-> Dirac::gamma_axis).";

ntGamma5::usage = "ntGamma5[din, dout] — Dirac gamma_5 on spinor axes din,dout.";

ntSigma::usage = "ntSigma[legA, legB, din, dout] — INTERNAL codegen token for the bare γ-commutator [A,B] on spinor axes din,dout, folded from the struct-7 quark-gluon-vertex σ^{μν}. Each leg is {\"slash\", mom} (a slashed momentum) or {\"free\", mu} (an open Lorentz id). Emitted by foldDiracSigma so the antisymmetric γ-pair is never distributed into two traces.";

ntDeltaDirac::usage = "ntDeltaDirac[din, dout] — spinor identity delta (Dirac::identity); closes a spinor trace.";

(* ---- frame builders (public): a frame fixes each momentum's four components as closed-form
        expressions in the runtime scalars, and is passed to NumTrace via "Frame"->frame. ---- *)

propFrame::usage = "propFrame[p, l1, cos1, q1, ql] — one external direction (q1 along axis 0) + loop (ql in the 0-1 plane at angle cos1). The Zc/ZA symmetric-point kinematics.";

sp3Frame::usage = "sp3Frame[p, l1, cos1, cos2, q1, q2, q3, ql] — three externals at the symmetric point (120 deg in the 0-1 plane) + loop with two angles. The A3 kinematics.";
gen3Frame::usage = "gen3Frame[p1m, p2m, cosP, l1, cos1, cos2, q1, q2, q3, ql] — three externals in a GENERAL configuration (independent magnitudes p1m,p2m and opening angle cosP; q3 = -(q1+q2)) + loop with two angles. Use instead of sp3Frame whenever the symmetric point would be a degenerate slice: the full AqbqDirect basis, for one, has Det[Gram] = 0 there. p1m==p2m, cosP==-1/2 recovers sp3Frame.";

sp4Frame::usage = "sp4Frame[p, l1, cos1, cos2, phi, q1, q2, q3, q4, ql] — four externals at a regular tetrahedron + loop with three angles. The A4 kinematics.";

propFrameFT::usage = "propFrameFT[p0, p, l0, l1, cos1, pSym, lSym] — finite-temperature frame: component 0 is the Matsubara/temporal direction (heat bath u=(1,0,0,0)). External pSym = {p0, p, 0, 0} (p0 e.g. the lowest fermionic Matsubara frequency pi T), loop lSym = {l0, l1 cos1, l1 sqrt(1-cos1^2), 0} with an independent temporal l0.";

sp3FrameFT::usage = "sp3FrameFT[{p10,p20,p30}, p, l0, l1, cos1, phi, q1, q2, q3, ql] — finite-temperature 3-point frame: three externals at the SPATIAL symmetric point (equal spatial length p, 120 deg apart in the spatial 1-2 plane, summing to zero), each carrying its own temporal component from p0s; loop ql = {l0, l1 sin1 cos(phi), l1 sin1 sin(phi), l1 cos1} with an independent Matsubara l0. Slot 0 is temporal, 1..3 spatial. The spatial configuration and the loop parametrisation match DiFfRG's DeclareSymmetricPoints3DP3 exactly, so cos(l,p_i) agrees term by term. The finite-T analogue of sp3Frame.";

sp4FrameFT::usage = "sp4FrameFT[{p10,p20,p30,p40}, p, l0, l1, cos1, phi, q1, q2, q3, q4, ql] — finite-temperature 4-point frame: four externals at the SPATIAL symmetric point (regular tetrahedron in the spatial slots, mutual cosines -1/3), each with its own temporal component; loop as in sp3FrameFT. Matches DiFfRG's DeclareSymmetricPoints3DP4. Note this needs only TWO loop angles where the vacuum sp4Frame needs three: four SPATIAL vectors summing to zero fit in three dimensions.";

frameShiftedLoop::usage = "frameShiftedLoop[frame, lSym, lfSym, shift] — add the FERMIONIC loop-momentum partner lfSym to a finite-T frame: same three-momentum as lSym, temporal component shifted by `shift` (typically pi T). FunKit routes the bare loop momentum onto the regulated line and names it after that line's statistics, so a quark self-energy carries BOTH l1 (gluon-regulator diagram) and lf1 (quark-regulator diagram) at once; with the integrator summing the bosonic tower f0 = 2 pi n T, lf0 = f0 + pi T is the odd tower and the shifted bosonic quadrature reproduces the fermionic sum exactly. Do NOT rewrite lf1 -> l1 at finite T (the vacuum derivations do, correctly, because there the two coincide).";

frameSpatialCosines::usage = "frameSpatialCosines[frame, ql, {q1,...,qn}, p, l1] — the loop-external SPATIAL cosines cos(l,p_i) read off `frame`, returned as the rule list {cosl1p1 -> ..., ...} that the dressing parametrisation and MakeNTKernel's \"AngleDefs\" both consume. Derive them rather than restating them: the tensor part resolves momenta through the frame while the dressings go through AngleDefs, and nothing downstream compares the two.";

Begin["`Private`"];

$NumTracerDirectory = DirectoryName[$InputFileName];

(* ---- generator thread / RAM caps ----
   The generator binary reads NT_GEN_MAXW (phase A: parallel Dirac-trace contraction) and
   NT_GEN_MAXW_B (phase B: per-net fold) via getenv at run time. SetEnvironment mutates the
   process environment that Run[]'s child shell inherits, so setting them here reaches the
   generator with no extra plumbing. Phase A's peak RSS is ~W x the per-contraction working
   set (a hard ceiling), and each phase-B worker holds one full contraction — so for dense
   flows nB is the dominant RAM lever. *)
SetNumTracerThreads[nA_Integer?Positive, nB_Integer?Positive] := (
   SetEnvironment["NT_GEN_MAXW" -> ToString[nA]];
   SetEnvironment["NT_GEN_MAXW_B" -> ToString[nB]];
   {nA, nB});

SetNumTracerThreads[n_Integer?Positive] := SetNumTracerThreads[n, n];

SetNumTracerThreads[] := GetNumTracerThreads[];

GetNumTracerThreads[] := Module[{a = Environment["NT_GEN_MAXW"], b = Environment["NT_GEN_MAXW_B"]},
   {If[a === $Failed, Automatic, ToExpression[a]],
    If[b === $Failed, Automatic, ToExpression[b]]}];


(* ---- load the implementation files, LOUDLY -------------------------------------------------
   A syntax error in one of these is the single most dangerous edit in this package, because
   Mathematica makes it a SILENT NO-OP: on `Syntax::sntue` Get returns Null (NOT $Failed — checking
   the return value does not work), every definition BEFORE the malformed expression is installed,
   and every definition after it — up to the end of the file, since an unbalanced bracket swallows
   the remainder into one unfinished expression — quietly is not.

   What that looks like downstream, measured: an unbalanced bracket in the middle of Codegen.m left
   `Options[MakeNTKernel]` undefined 1800 lines later, so a generator script reported
   `SetOptions::optnf: RuntimeInclude is not a known option for MakeNTKernel`, ran its numeric
   backend in 0.0001 s, emitted NOTHING, printed "kernels generated" and exited 0. The real cause
   was one line in the load output, thousands of lines earlier, and nothing connected the two.

   Check with an explicit message list turns that into an abort at load time. It is inert on a clean
   file (verified against both an unfinished expression and a stray closing bracket). *)
ntLoadPart[file_String] :=
  If[Check[Get[FileNameJoin[{$NumTracerDirectory, file}]], $Failed,
       Syntax::sntue, Syntax::sntx, Syntax::sntxi,
       Syntax::bktmcp, Syntax::bktmop, Syntax::tsntxi] === $Failed,
    Message[NumTrace::loadsyntax, file]; Abort[]];

ntLoadPart["DSL.m"];

ntLoadPart["Frames.m"];

ntLoadPart["Codegen.m"];

ntLoadPart["FunKitAdapter.m"];

ntLoadPart["DiFfRG_compat.m"];

End[];

Protect["NumTracer`*"];

EndPackage[];
