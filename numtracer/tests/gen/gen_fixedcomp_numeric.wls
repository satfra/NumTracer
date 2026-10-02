#!/usr/bin/env wolframscript
(* Codegen gate for FIXED Lorentz components (gamma^0, gamma^1, ... — the finite-T 3+1 split in
   which the four-quark Fierz bases are written).

   A fixed component is a Lorentz index pinned to a CONCRETE value: ntGamma[0, d1, d2] is gamma^0,
   NOT "gamma^mu with a label named 0". Before the fix, DSL.m's label machinery had no such notion
   (labelsOf[ntGamma[mu,...]] returns {mu,...} unconditionally), so the literal 0 was read as an
   ordinary contraction label, occurred in many factors at once, and checkLabels aborted:

       NumTrace::privclash: label(s) 0 are private dummies of one factor
                            but also occur outside it

   That guard was behaving CORRECTLY — the feature was simply missing. It is now supplied by a
   rewrite in NumTrace (expandFixedComponents): a fixed component is a contraction with the constant
   unit basis vector e_i, so gamma^i becomes an ordinary SLASH against e_i and the whole thing needs
   no new DFac kind and no C++ change at all.

   The frame here is propFrameFT, whose TEMPORAL components are nonzero — at a T=0 vacuum frame
   every q_0 vanishes and each oracle below would collapse to 0, i.e. the gate would pass vacuously
   no matter what the tracer did. *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3;

args = {Global`p0, Global`p, Global`l0, Global`l1, Global`cos1};

(* q1 = {p0, p, 0, 0};  ql = {l0, l1 cos1, l1 sin1, 0} *)
frame = NumTracer`propFrameFT[Global`p0, Global`p, Global`l0, Global`l1, Global`cos1,
                              Global`q1, Global`ql];

(* The REFERENCE frame adds an explicit unit vector as an ordinary user momentum, so the same trace
   can be written WITHOUT the fixed-component feature: gamma^0 == gamma^mu contracted with e0_mu.
   Grading the two against each other is convention-independent — it cannot be fooled by a sign or
   metric-convention error in the hand-written oracles below, which is the failure mode that has
   twice produced a false green in this codebase. *)
frameRef = Append[frame, Global`e0 -> {1, 0, 0, 0}];

emit[net_, nm_, ns_, fr_] := Module[{ntk},
  ntk = NumTrace[net, "Frame" -> fr, "Args" -> args, "Dressings" -> {}];
  If[Length[ntk[[1]]["Diagrams"]] < 1,
    Print["FAIL: net ", nm, " produced 0 diagrams (was it pruned as an odd trace?)"]; Exit[1]];
  MakeNTKernel[ntk,
    FileNameJoin[{$gen, "gen_" <> ns <> ".cpp"}],
    FileNameJoin[{$gen, nm <> ".hh"}],
    FileNameJoin[{$gen, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> {}];
  Print["[gen] ok: ", nm]];

(* --- Assertion A: the rewrite happens at all, and produces a SLASH against a unit vector. -------
   Pre-fix this NumTrace call aborts outright with privclash. Post-fix the fixed index must be gone
   from the label machinery and the unit vector must have reached the env as a momentum (without an
   env Base, compileDirac's slash emission cannot reference it). *)
netA = ntGamma[0, d1, d2] ntGamma[a1, d2, d1] ntVec[Global`q1, a1];
ntkA = NumTrace[netA, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[! MemberQ[Keys[ntkA[[1]]["Env"]], ntUnitVec[0]],
  Print["FAIL: the unit basis vector e_0 did not reach the env; a fixed-component gamma cannot ",
        "be emitted as a slash without one. Env = ", Keys[ntkA[[1]]["Env"]]]; Exit[1]];
Print["[gen] ok: fixed component rewritten, e_0 in env"];

(* --- Assertion B: 0-based component indexing. -----------------------------------------------
   tr(gamma^0 q1_slash) = 4 q1_0 = 4 p0   and   tr(gamma^1 q1_slash) = 4 q1_1 = 4 p.
   In propFrameFT p0 and p are INDEPENDENT args, so an off-by-one in the unit vector (e_i placed at
   slot i instead of i+1, or 1-based components) swaps the two answers and is caught. This is the
   whole reason both are emitted rather than just one. *)
emit[ntGamma[0, d1, d2] ntGamma[a1, d2, d1] ntVec[Global`q1, a1], "Fixc_t_kernel", "fixc_t", frame];
emit[ntGamma[1, d1, d2] ntGamma[a1, d2, d1] ntVec[Global`q1, a1], "Fixc_x_kernel", "fixc_x", frame];

(* --- Assertion C: a fixed component contracted with ITSELF. ----------------------------------
   tr(gamma^0 gamma^0) = 4 (e_0 . e_0) = 4, a pure constant. Two fixed components on the same chain
   get two DISTINCT fresh dummies but the SAME env momentum, so this pins that the rewrite neither
   fuses the two dummies (which would make the pair a self-contraction) nor allocates e_0 twice. *)
emit[ntGamma[0, d1, d2] ntGamma[0, d2, d1], "Fixc_tt_kernel", "fixc_tt", frame];

(* --- Assertion D: the full four-gamma trace, fixed components INTERLEAVED with ordinary slashes.
   tr(gamma^0 q1_slash gamma^0 ql_slash) = 4 [2 q1_0 ql_0 - (q1.ql)] = 4 [p0 l0 - p l1 cos1].
   Kinematics-dependent in every argument, so it cannot fold to a constant. *)
netD = ntGamma[0, d1, d2] ntGamma[a1, d2, d3] ntVec[Global`q1, a1] *
       ntGamma[0, d3, d4] ntGamma[a2, d4, d1] ntVec[Global`ql, a2];
emit[netD, "Fixc_4_kernel", "fixc_4", frame];

(* --- Assertion E: the convention-independent cross-check. ------------------------------------
   The SAME trace as D, but with gamma^0 written the long way as gamma^mu contracted against an
   explicit unit-vector momentum supplied by the user's frame. This uses none of the new code path,
   so agreement between D and E proves the rewrite is faithful independently of whether my oracle's
   Euclidean sign conventions are right. *)
netE = ntGamma[b1, d1, d2] ntVec[Global`e0, b1] ntGamma[a1, d2, d3] ntVec[Global`q1, a1] *
       ntGamma[b2, d3, d4] ntVec[Global`e0, b2] ntGamma[a2, d4, d1] ntVec[Global`ql, a2];
emit[netE, "Fixc_4ref_kernel", "fixc_4ref", frameRef];

(* --- Assertion F: g^{0 nu} — a fixed component on a METRIC, not a gamma. ----------------------
   g^{0 nu} q1_nu = q1_0 = p0, so tr(gamma^0 ql_slash) g^{0 nu} q1_nu = 4 ql_0 p0 = 4 l0 p0. *)
netF = ntMetric[0, nu] ntVec[Global`q1, nu] *
       ntGamma[0, d1, d2] ntGamma[a1, d2, d1] ntVec[Global`ql, a1];
emit[netF, "Fixc_g_kernel", "fixc_g", frame];

(* --- Assertion G: a fixed component on a SIGMA leg. -------------------------------------------
   sigma's free leg carries an OPEN Lorentz id; a fixed component is not a free leg, so the rewrite
   must convert it to the existing SLASH leg kind instead. tr(sigma^{0 nu} ...) exercises the
   comm overloads with a slashed leg rather than two free legs. *)
netG = ntSigma[{"free", 0}, {"slash", {{1, Global`q1}}}, d1, d2] *
       ntGamma[a1, d2, d3] ntVec[Global`ql, a1] *
       ntGamma[a2, d3, d1] ntVec[Global`q1, a2];
ntkG = NumTrace[netG, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[! FreeQ[ntkG[[1]]["Diagrams"], {"free", _Integer}],
  Print["FAIL: a sigma leg still carries a fixed component as a FREE leg; it must be rewritten ",
        "to a slash leg against the unit basis vector."]; Exit[1]];
Print["[gen] ok: sigma fixed-component leg rewritten to a slash leg"];

Print["Fixed-component (gamma^0) codegen kernels generated."];
