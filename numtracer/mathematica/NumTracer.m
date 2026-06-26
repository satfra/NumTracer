(* ::Package:: *)

(* NumTracer` — Mathematica front-end that turns a traced tensor network into a
   C++ kernel driving the compile-time `et` engine (include/numtracer/).

   Two public entry points:
     NumTrace[net, opts]     builds an intermediate-expression tree (an NTKernel)
                             from a network written in the DSL heads below;
     MakeNTKernel[ntk, ...]  serialises that tree to C++; "Backend" selects the
                             path ("ET" | "Dense" -> one kernel header, "Generate"
                             -> the invariant-basis build-time generation path).

   The network is `Σ_terms coeff(dressings…) × trace(momenta)`: the trace is a
   contraction of the tensor heads (emitted as `et` builder calls + contract_all,
   lowered by expr::eval_real), the coeff is flat C++. This is exactly the split
   tests/refshim/flow_ym_zc.hpp demonstrates by hand — we generate it. *)

BeginPackage["NumTracer`"];

Unprotect["NumTracer`*"];

ClearAll["NumTracer`*"];

(* ---- public entry points ---- *)

NumTrace::usage = "NumTrace[net, \"Frame\"->frame, \"Args\"->{...}] analyses a DSL tensor network into an NTKernel[<|...|>]: a list of terms, each a scalar coefficient times one or more independent contraction components, plus the env-id layout and frame needed to emit code.";

MakeNTKernel::usage = "MakeNTKernel[ntk, file, \"Backend\"->\"ET\"|\"Dense\", \"Name\"->\"X\", \"Dressings\"->{...}] serialises an NTKernel to a single C++ kernel header at `file` (symbolic ETensor kernel, or the brute-force numeric DTensor baseline).\nMakeNTKernel[ntk, genFile, kernelFile, tracesFile, \"Backend\"->\"Generate\", \"Name\"->\"X\", \"Namespace\"->ns, ...] emits the invariant-basis C++ generation path: a build-time generator program at `genFile`, runs it to produce the committed straight-line traces header `tracesFile`, and the kernel header `kernelFile` that fills the fundamental symbols and calls them.\n\"Backend\"->Automatic (default) resolves from the argument count: one output file -> \"ET\", three -> \"Generate\".";

NTKernel::usage = "NTKernel[assoc] is the analysed intermediate-expression tree produced by NumTrace and consumed by MakeNTKernel.";

FromFunKit::usage = "FromFunKit[flow, \"FlavourGroup\"->n] rewrites a FunKit traced flow (after // dressingRules) into the NumTracer DSL; scalar products become ntSP, resolved by the frame. SU(N) group tokens get their rank baked into the ntSUN* heads — colour from Global`Nc, the isospin group from \"FlavourGroup\" (default Global`Nf or 2).";

MakeNTKernelDiFfRG::usage = "MakeNTKernelDiFfRG[ntk, \"Name\"->\"ZA\", \"Integrator\"->\"Integrator_p2_1ang\", \"Parameters\"->kernelParameterList, \"IntegrationVariables\"->{\"l1\",\"cos1\"}, \"AngleDefs\"->{...}] scaffolds the DiFfRG flows/<Name>/ plumbing (MakeKernel on a placeholder body) and overwrites kernel.hh/kernels.hh with the NumTracer-traced numeric kernel. Dressings and their interpolator type are auto-derived from Parameters; the Regulator alias and the DiFfRG emission constants (namespace, includes, support API) are baked in. A loop-independent term flat-added to the integral (DiFfRG's constExpr) can be passed either positionally — MakeNTKernelDiFfRG[ntk, constExpr, opts] — or as \"Constant\"->expr; it populates constant(p,k,dressings) and is a plain expression (e.g. ZA[p]), not an NTKernel. Pair with UpdateNTFlows. Requires DiFfRG to be loaded.";

UpdateNTFlows::usage = "UpdateNTFlows[name] runs DiFfRG's UpdateFlows[name] then idempotently patches flows/CMakeLists.txt (find_package(NumTracer) + NumTracer::NumTracer/constexpr_budget link + UNITY_BUILD OFF). Bundling the two makes the CMake patch atomic, so it can never be left un-applied after an UpdateFlows regenerates the file. Requires DiFfRG to be loaded.";

(* ---- DSL tensor heads (inert tags; first arg of the momentum-bearing ones is
        a momentum/4-vector symbol, the rest are contraction index labels) ---- *)

ntMetric::usage = "ntMetric[mu, nu] — Euclidean Lorentz metric delta_{mu nu}.";

ntVec::usage = "ntVec[q, mu] — tensor leg q_mu of momentum q (mu a symbolic Lorentz label). ntVec[q, i] with a literal INTEGER i (0-based, 0 = temporal/Matsubara) is the scalar component q_i, resolved by the frame like ntSP — used for finite-T objects such as ntVec[p, 0] = pi T.";

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

ntSUNDiagFund::usage = "ntSUNDiagFund[N, i, j, name, scale] — a fundamental Kronecker delta_{ij} (rank N) carrying a PER-COMPONENT dressing: the diagonal insertion diag(name[0],..,name[N-1]) with each component evaluated at the kinematic `scale`. `name` is a kernel parameter (a std::array of N interpolators). Folds (via sun_value_dressed) to Σ_i name[i](scale) instead of the flavour-blind δ — this is how the u- and d-quark (etc.) are dressed differently WITHIN the SU(N) trace, no per-flavour diagram split.";

ntSUNDiagAdj::usage = "ntSUNDiagAdj[N, a, b, name, scale] — an adjoint Kronecker delta^{ab} (rank N) carrying a PER-COMPONENT adjoint dressing diag(name[0],..,name[N^2-2]), each evaluated at `scale`. `name` is a kernel parameter (a std::array of N^2-1 interpolators). Folds to Σ_a c_a name[a](scale) — e.g. a condensed gluon with a colour-component-dependent dressing Z_A^a.";

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

sp4Frame::usage = "sp4Frame[p, l1, cos1, cos2, phi, q1, q2, q3, q4, ql] — four externals at a regular tetrahedron + loop with three angles. The A4 kinematics.";

propFrameFT::usage = "propFrameFT[p0, p, l0, l1, cos1, pSym, lSym] — finite-temperature frame: component 0 is the Matsubara/temporal direction (heat bath u=(1,0,0,0)). External pSym = {p0, p, 0, 0} (p0 e.g. the lowest fermionic Matsubara frequency pi T), loop lSym = {l0, l1 cos1, l1 sqrt(1-cos1^2), 0} with an independent temporal l0.";

Begin["`Private`"];

$NumTracerDirectory = DirectoryName[$InputFileName];

Get[FileNameJoin[{$NumTracerDirectory, "DSL.m"}]];

Get[FileNameJoin[{$NumTracerDirectory, "Frames.m"}]];

Get[FileNameJoin[{$NumTracerDirectory, "Codegen.m"}]];

Get[FileNameJoin[{$NumTracerDirectory, "FunKitAdapter.m"}]];

Get[FileNameJoin[{$NumTracerDirectory, "DiFfRG_compat.m"}]];

End[];

Protect["NumTracer`*"];

EndPackage[];
