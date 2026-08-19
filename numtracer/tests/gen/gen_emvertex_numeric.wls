#!/usr/bin/env wolframscript
(* Codegen gate for the FINITE-T ELECTRIC/MAGNETIC VERTEX SPLIT — the front-end path that
   arXiv:1708.03482 App. B/D puts under the three- and four-gluon vertices, and that
   qcd-codes/finite_T/YangMills/SP_EM now depends on.

   WHAT THIS COVERS THAT gen_ftproj_numeric.wls DOES NOT. That gate exercises one electric and one
   magnetic projector, both on the SAME loop momentum, each sandwiched between two vectors. The
   vertex split needs something structurally harder: THREE projectors on THREE DIFFERENT external
   momenta, contracted through a shared multi-index core, and summed over which legs are electric.
   That is where an index/momentum mix-up (the projector on leg i must carry leg i's own momentum)
   or a lost Pi^M Pi^E = 0 stops being visible.

   THE IDENTITIES. Write V[s] for the classical three-gluon Lorentz core with projector s_i attached
   to leg i, s_i in {M, E, T}. Two exact statements carry the whole (D1) projection argument:

     (a) ORTHOGONALITY   <V[MMM], V[EEM+perm]> = 0
         — leg by leg, because Pi^M(q) Pi^E(q) = 0. Without it the TensorBases dual projector is not
           the paper's T_i/<T_i,T_i> and every dressing extracted downstream is a mixture.
     (b) COMPLETENESS    <V[TTT], V[TTT]> = sum over all 8 assignments of <V[s], V[s]>
         — because Pi^M + Pi^E = Pi_perp per leg AND all the cross terms vanish by (a). This is the
           identity that lets the two-term compact vertex reproduce the pre-split one exactly.

   GRADING. Both are ZERO statements, and a zero kernel passes a zero test for free — the vacuous
   pass this repo has been bitten by before. So a third, deliberately NON-zero kernel is emitted
   alongside (the diagonal <V[MMM],V[MMM]> + <V[EEM],V[EEM]>), and compare_emvertex_num.cpp requires
   it to be large before it believes the zeros. The zeros themselves are graded ABSOLUTELY against
   that scale, never pointwise-relative — two roundings of zero always disagree by 100%.

   WHERE THE INDEPENDENCE COMES FROM. There is no dense-oracle twin here: MakeNTKernel no longer
   accepts "Backend" -> "Dense" — and now REFUSES an unknown option name rather than swallowing it,
   which is what let three stale generators go on passing it (they have since been deleted, along
   with the orphaned dense kernels they emitted). Identity (b) supplies the
   independence instead, and more sharply than dense-vs-numeric would: its left-hand side is built
   from the TRANSVERSE leaf, which carries 1/q^2 atoms, and its right-hand side from the electric
   and magnetic leaves, which carry 1/|q_vec|^2. Two structurally different lowerings of the same
   number, required to agree exactly. A bug that moved both in step would have to respect that.

   FRAME. sp3FrameFT with NONZERO external Matsubara components. At p0 = 0 the electric projector
   degenerates to the rank-1 u (x) u and half of these statements become true for the wrong reason;
   the same lesson as gen_ortho_numeric.wls, which had to leave the symmetric point because the
   basis is degenerate exactly there. *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

(* Three external legs at the spatial symmetric point, each with its OWN Matsubara component, plus
   the loop. p0s is deliberately {p01, p02, -p01-p02}-like: momentum conservation is not imposed on
   the temporal slots by the frame, and these identities hold leg-wise regardless, so three
   independent frequencies give the projectors the most generic kinematics available here. *)
frame = NumTracer`sp3FrameFT[{Global`p01, Global`p02, Global`p03}, Global`p,
          Global`f0, Global`l1, Global`cos1, Global`phi,
          Global`q1, Global`q2, Global`q3, Global`ql];
args  = {Global`p01, Global`p02, Global`p03, Global`p, Global`f0, Global`l1, Global`cos1, Global`phi};

(* projector leaf by label; each carries its own leg's momentum *)
pr["M", q_, i_, j_] := ntMagneticProj[q, i, j];
pr["E", q_, i_, j_] := ntElectricProj[q, i, j];
pr["T", q_, i_, j_] := ntTransProj[q, i, j];

legs = {Global`q1, Global`q2, Global`q3};

(* The classical three-gluon Lorentz core, indices supplied. Colour is dropped on purpose: the E/M
   split is a Lorentz statement, and an SU(N) factor would only make the gate slower. The loop
   momentum ql rides in through the momentum differences so the kernel genuinely depends on the
   loop angles rather than collapsing to a function of p alone. *)
core[m_, n_, r_] :=
  ntMetric[m, n] ntVec[Global`q2 - Global`q1 + Global`ql, r] +
  ntMetric[m, r] ntVec[Global`q1 - Global`q3 - Global`ql, n] +
  ntMetric[n, r] ntVec[Global`q3 - Global`q2, m];

(* V[s] with a caller-chosen index triple, so two copies can be contracted without collision *)
V[sel_List, {m_, n_, r_}] := Module[{ii = {m, n, r}, jj},
  jj = Table[Unique["vi"], {3}];
  Times @@ MapThread[pr[#1, #2, #3, #4] &, {sel, legs, ii, jj}] core @@ jj];

(* <V[s], V[s']> : contract the three outer indices *)
ip[s1_List, s2_List] := Module[{o = Table[Unique["oi"], {3}]},
  V[s1, o] V[s2, o]];

symV[sel_List] := DeleteDuplicates[Permutations[sel]];
ipSym[s1_List, s2_List] := Total[Outer[ip[#1, #2] &, symV[s1], symV[s2], 1], 2];

allAssignments = Tuples[{"M", "E"}, 3];

nets = <|
  (* NON-ZERO control, and the numeric-vs-dense engine check *)
  "diag" -> ipSym[{"M", "M", "M"}, {"M", "M", "M"}] + ipSym[{"E", "E", "M"}, {"E", "E", "M"}],
  (* (a) orthogonality — must vanish identically *)
  "off"  -> ipSym[{"M", "M", "M"}, {"E", "E", "M"}],
  (* (b) completeness — must vanish identically *)
  "comp" -> ip[{"T", "T", "T"}, {"T", "T", "T"}] -
            Total[ip[#, #] & /@ allAssignments]
|>;

emit[key_String, nm_String, ns_String] := Module[{ntk},
  ntk = NumTrace[nets[key], "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
  Print["[emvertex] ", key, ": ", Length[ntk[[1]]["Diagrams"]], " diagram(s)"];
  (* A net that traced to NOTHING would make every downstream assertion vacuously true. The two
     zero-identities are *expected* to evaluate to zero, but they must do so by cancellation of real
     diagrams, not by producing none. *)
  If[Length[ntk[[1]]["Diagrams"]] < 1,
    Print["FAIL: net \"", key, "\" produced NO diagrams — the identity would then be confirmed by ",
          "an empty kernel rather than by cancellation."]; Exit[1]];
  MakeNTKernel[ntk,
    FileNameJoin[{$gen, "gen_emvertex_" <> ns <> ".cpp"}],
    FileNameJoin[{$gen, nm <> ".hh"}],
    FileNameJoin[{$gen, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> {}];
  ];

emit["diag", "EMVert_diag_kernel", "emvert_diag"];
emit["off",  "EMVert_off_kernel",  "emvert_off"];
emit["comp", "EMVert_comp_kernel", "emvert_comp"];

(* The E/M leaves must actually have reached the emitted generator source. Without this the gate
   could pass on a net where the projectors were silently rewritten away. *)
gsrc = StringJoin[(Import[#, "Text"] &) /@
   Select[FileNames["gen_emvertex_*.cpp", $gen], FileExistsQ]];
If[! (StringContainsQ[gsrc, "eproj<"] && StringContainsQ[gsrc, "mproj<"]),
  Print["FAIL: emitted generator source missing eproj<...>() / mproj<...>() leaves"]; Exit[1]];
Print["[emvertex] ok: eproj/mproj leaves present in the generator source."];
Print["E/M vertex-split codegen kernels generated."];
