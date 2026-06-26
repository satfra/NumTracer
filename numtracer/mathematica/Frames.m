(* ::Package:: *)

(* Frame handling: a frame fixes each momentum's four components as closed-form
   expressions in the runtime scalars (l1, cos1, cos2, p). The trace is then a
   network over Var leaves filled per grid point; the frame is the generated
   analogue of lorentz.hpp's along0 / at_angle.

   A frame is an association  <| qvec -> {c0, c1, c2, c3}, ... |>  keyed by the
   (vector) momentum symbols that appear in ntVec/ntTransProj, with values in the
   scalar args. Composite momenta resolve by list arithmetic, so a difference like
   q1 - ql is handled automatically. *)

(* Four components of a (possibly composite) momentum in the frame. *)

resolveComponents[q_, frame_] :=
  Simplify[q /. Normal[frame]];

(* Euclidean q^2 of a momentum in the frame (drives the 1/q^2 env slot). *)

normSqExpr[q_, frame_] :=
  Simplify[Total[resolveComponents[q, frame] ^ 2]];

(* Spatial |q⃗|^2 of a momentum in the frame: drop the temporal component (slot 0) and sum the
   squares of the spatial three. Drives the 1/|q⃗|^2 env slot of the finite-T electric/magnetic
   projectors. At T=0 vacuum frames (all temporal slots zero) this coincides with normSqExpr. *)

spatialNormSqExpr[q_, frame_] :=
  Simplify[Total[Rest[resolveComponents[q, frame]] ^ 2]];

(* ---- symmetric-point frames (reproduce the FormTracer SP kernels' kinematics) ----

   The external gluons sit at 120 degrees in the 0-1 plane (|p_i| = p, sum = 0); the
   loop l makes angle cos1 with q1 and is parametrised by a second angle cos2. These
   reproduce the oracle's cosl1p2/p3 = (-cos1 +- sqrt(3-3 cos1^2) cos2)/2 exactly.

   Builders take the scalar arg symbols explicitly to stay context-clean; the caller
   chooses the vector tags (here q1,q2,q3,ql) and uses the same in its DSL network. *)

(* One external direction + loop: q1 along axis 0, l in the 0-1 plane. (Zc / ZA) *)

propFrame[p_, l1_, cos1_, q1_, ql_] :=
  <|q1 -> p {1, 0, 0, 0}, ql -> l1 {cos1, Sqrt[1 - cos1^2], 0, 0}|>;

(* Three externals at the symmetric point + loop with two angles. (A3) *)

sp3Frame[p_, l1_, cos1_, cos2_, q1_, q2_, q3_, ql_] :=
  <|q1 -> p {1, 0, 0, 0}, q2 -> p {-1/2, Sqrt[3] / 2, 0, 0}, q3 -> p {-1/2, -Sqrt[3] / 2, 0, 0}, ql -> l1 {cos1, Sqrt[1 - cos1^2] cos2, Sqrt[(1 - cos1^2) (1 - cos2^2)], 0}|>;

(* Four externals at the symmetric point + loop with three angles. (A4)

   The four external gluons sit at the vertices of a regular tetrahedron in the 3D
   subspace (slots 0-2): |p_i| = p, sum = 0, p_i.p_j = -p^2/3. The loop l makes angle
   cos1 with the perpendicular (slot-3) axis; its in-plane part (magnitude
   sqrt(1-cos1^2)) is parametrised by cos2 (toward q1) and phi (rotation in the plane
   perpendicular to q1, spanned by the orthonormal perpA,perpB below). This reproduces
   the oracle's DeclareSymmetricPoints4DP4 loop-external cosines cosl1p{1..4} exactly
   (verified against SP4Defs). *)

sp4Frame[p_, l1_, cos1_, cos2_, phi_, q1_, q2_, q3_, q4_, ql_] :=
  Module[{v1, v2, v3, v4, perpA, perpB, sin1, sin2, loopDir},
    v1 = {1, 1, 1} / Sqrt[3];
    v2 = {1, -1, -1} / Sqrt[3];
    v3 = {-1, 1, -1} / Sqrt[3];
    v4 = {-1, -1, 1} / Sqrt[3];(* tetrahedron vertices *)
    perpA = {2, -1, -1} / Sqrt[6];
    perpB = {0, 1, -1} / Sqrt[2];(* orthonormal, perpendicular to v1 *)
    sin1 = Sqrt[1 - cos1^2];
    sin2 = Sqrt[1 - cos2^2];
    loopDir = cos2 v1 + sin2 (Cos[phi] perpA + Sin[phi] perpB); (* loop in-plane direction *)
    <|q1 -> p Append[v1, 0], q2 -> p Append[v2, 0], q3 -> p Append[v3, 0], q4 -> p Append[v4, 0], ql -> l1 Append[sin1 loopDir, cos1]|>
  ];

(* ---- finite-temperature frame (O(4) -> O(3) broken) -------------------------------------------

   Component 0 is the Euclidean-time / Matsubara direction (the heat-bath rest frame u=(1,0,0,0));
   1..3 are spatial. The external momentum carries an independent temporal component p0 (e.g. the
   lowest fermionic Matsubara frequency πT) and a spatial magnitude p along axis 1; the loop carries
   an independent temporal component l0 (the summed Matsubara frequency) and a spatial 3-vector of
   magnitude l1 making angle cos1 with the external spatial momentum (azimuthal symmetry puts it in
   the 1-2 plane, slot 3 unused). This is the finite-T analogue of propFrame, with the temporal
   components no longer zero — so sps/vec[·,0] and the electric/magnetic projectors are meaningful. *)

propFrameFT[p0_, p_, l0_, l1_, cos1_, pSym_, lSym_] :=
  <|pSym -> {p0, p, 0, 0}, lSym -> {l0, l1 cos1, l1 Sqrt[1 - cos1^2], 0}|>;
