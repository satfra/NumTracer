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

resolveComponents::zeroMomentum = "A momentum in the frame resolved to `1` instead of a list of four \
components. This means a tensor leaf carries the ZERO four-vector -- almost always because a \
kinematic substitution (p -> 0 for a p^2-coefficient extraction, a soft-meson slice, ...) was \
applied to an expression that still contained a momentum-dependent PROJECTOR. transProj[0,mu,nu] is \
0/0, so there is nothing sensible to resolve. Keep the projector out of the substitution (subtract \
inside the diagram only) or pick a frame in which the leg is nonzero.";

resolveComponents[q_, frame_] :=
  Module[{c = Simplify[q /. Normal[frame]]},
(* Guard rather than let a bare 0 travel on: downstream this becomes `mpcpp /@ 0`, which returns the
   INTEGER 0 where a C++ string was expected, and the whole generator source degenerates into one
   unevaluated StringJoin. The symptom is then a StringTake::strse plus a several-hundred-line dump
   of the half-built source -- true, but useless. Fail here, where the cause is still visible. *)
    If[! (ListQ[c] && Length[c] === 4),
      Message[resolveComponents::zeroMomentum, c];
      Abort[]];
    c];

(* `normSqExpr` (Euclidean q^2) and `spatialNormSqExpr` (spatial |q⃗|^2) used to live here, to build
   the 1/q^2 and 1/|q⃗|^2 env slots. Both are now derived where the env is built, from the same
   `resolveComponents` output the component table already holds, so neither had a caller. *)

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

(* Three externals in a GENERAL configuration + loop with two angles.

   sp3Frame pins the symmetric point (|p_i| all equal, p_i.p_j = -p^2/2), which is a MEASURE-ZERO
   slice of the 3-point phase space — and some tensor bases are degenerate exactly there. The full
   AqbqDirect basis is the case in point: its 12 structures are linearly DEPENDENT at the symmetric
   point (Det[TBGetMetric["AqbqDirect"]] = 0), so its inverse metric is 0/0 and every dual projector
   built from it comes out Indeterminate. A basis-orthonormality check must therefore be traced in a
   GENERAL frame, or it is only testing the tracer on a special point where the algebra degenerates.

   Two independent external magnitudes p1m, p2m and the angle cosP between them; the third leg is
   fixed by momentum conservation q3 = -(q1+q2), so |q3| varies independently of |q1|,|q2|. The loop
   keeps sp3Frame's parametrisation (cos1 to axis 0, cos2 in the 0-1 plane). Setting
   p1m == p2m and cosP == -1/2 recovers the symmetric point exactly. *)

gen3Frame[p1m_, p2m_, cosP_, l1_, cos1_, cos2_, q1_, q2_, q3_, ql_] :=
  Module[{v1, v2},
    v1 = p1m {1, 0, 0, 0};
    v2 = p2m {cosP, Sqrt[1 - cosP^2], 0, 0};
    <|q1 -> v1, q2 -> v2, q3 -> -(v1 + v2),
      ql -> l1 {cos1, Sqrt[1 - cos1^2] cos2, Sqrt[(1 - cos1^2) (1 - cos2^2)], 0}|>];

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

(* ---- finite-temperature SPATIAL symmetric-point frames ---------------------------------------

   The finite-T analogues of sp3Frame / sp4Frame. Once O(4) breaks to O(3) the "symmetric point"
   is a statement about the SPATIAL momenta only: the n external three-momenta have equal length p,
   sum to zero, and sit at mutual angles arccos(-1/(n-1)). Each leg additionally carries its own
   independent temporal (Matsubara) component, passed as the list `p0s` -- typically all zero for
   gluon legs, but the list form also covers legs pinned to a nonzero external frequency.

   The spatial configurations reproduce DiFfRG's DeclareSymmetricPoints3DP3 / ...3DP4 EXACTLY,
   including the loop parametrisation: with vec3[T,F] = {Sin[T] Cos[F], Sin[T] Sin[F], Cos[T]} the
   loop direction is vec3[ArcCos[cos1], phi], so cos1 is the polar angle against the 3-axis and phi
   the azimuth. Consequently the loop-external cosines cos(l,p_i) = loopDir . V_i agree with
   DeclareSymmetricPoints3D* term by term, and a kernel generated here can be checked against a
   FormTracer kernel that used those declarations.

   Note the loop's spatial part is a UNIT 3-vector times l1 -- Sum_{i=1..3} U_i^2 = 1 -- which is
   precisely the property unitLoopSpatialQ / unitLoopMixedFrameSpec in Codegen.m need to collapse
   the bare-loop denominator to the two-term l0^2 + l1^2. Do not reparametrise the loop in a way
   that breaks that identity (e.g. by inlining a non-polar direction), or the emitted traces blow
   up. Slot 0 is temporal throughout, slots 1..3 spatial. *)

(* Three externals at the spatial symmetric point (120 degrees in the spatial 1-2 plane) + loop
   with polar angle cos1 (against the spatial 3-axis) and azimuth phi. (ZAcbc / ZA3) *)

sp3FrameFT[p0s_List, p_, l0_, l1_, cos1_, phi_, q1_, q2_, q3_, ql_] :=
  <|q1 -> {p0s[[1]], p, 0, 0},
    q2 -> {p0s[[2]], -p/2, p Sqrt[3]/2, 0},
    q3 -> {p0s[[3]], -p/2, -p Sqrt[3]/2, 0},
    ql -> {l0, l1 Sqrt[1 - cos1^2] Cos[phi], l1 Sqrt[1 - cos1^2] Sin[phi], l1 cos1}|>;

(* Four externals at the spatial symmetric point (regular tetrahedron in the spatial slots, mutual
   cosines -1/3) + the same loop parametrisation. (ZA4)

   Unlike the vacuum sp4Frame -- where four momenta summing to zero need all four Euclidean slots
   and hence a THIRD loop angle -- four SPATIAL vectors summing to zero fit in three dimensions, so
   the finite-T ZA4 kernel needs only (cos1, phi) plus the Matsubara variable. That is why the
   finite-T ZA4 uses Integrator_fT_p2_4D_2ang, not a 3-angle integrator. *)

sp4FrameFT[p0s_List, p_, l0_, l1_, cos1_, phi_, q1_, q2_, q3_, q4_, ql_] :=
  Module[{s = 2 Sqrt[2]/3},(* Sin[ArcCos[-1/3]] *)
    <|q1 -> {p0s[[1]], 0, 0, p},
      q2 -> {p0s[[2]], p s, 0, -p/3},
      q3 -> {p0s[[3]], -p s/2, p s Sqrt[3]/2, -p/3},
      q4 -> {p0s[[4]], -p s/2, -p s Sqrt[3]/2, -p/3},
      ql -> {l0, l1 Sqrt[1 - cos1^2] Cos[phi], l1 Sqrt[1 - cos1^2] Sin[phi], l1 cos1}|>];

(* ---- the fermionic loop-momentum partner (finite T) ------------------------------------------

   At T > 0 a diagram can carry TWO loop-momentum tags at once. FunKit routes the BARE loop
   momentum onto whichever line carries the regulator insertion, and names it after that line's
   statistics: `l1` when d_t R sits on a Bose line, `lf1` when it sits on a Fermi line (see
   FunKit modules/AnSEL.m). A quark self-energy therefore produces both — the gluon-regulator
   diagram routed through l1 and the quark-regulator diagram through lf1 — inside one expression.

   The two are the SAME three-momentum; they differ only in which Matsubara tower the temporal
   component runs over. With the integrator summing the BOSONIC tower f0 = 2 pi n T, the fermionic
   partner is that same variable shifted by pi T, so that l0 stays even and lf0 = f0 + pi T is odd.
   Shifting the bosonic quadrature by pi T reproduces the fermionic sum exactly:
   T sum_n g(2 pi n T) with g(x) = h(x + pi T) IS T sum_n h((2n+1) pi T).

   Adding the partner as a separate frame key (rather than rewriting lf1 -> l1, which is what the
   vacuum derivations do and what would be WRONG here) keeps the two towers distinct while the
   spatial parts stay literally identical, so composite momenta like lf1 - p1 resolve correctly by
   the frame's ordinary list arithmetic.

   `shift` is added to component 0 only; passing -pi T gives the other sign convention. *)

frameShiftedLoop[frame_, lSym_, lfSym_, shift_] :=
  Join[frame, <|lfSym -> MapAt[# + shift &, frame[lSym], 1]|>];

(* Loop-external SPATIAL cosines cos(l, p_i) read straight off a frame, as the rule list the
   dressing parametrisation (SPParam) and MakeNTKernel's "AngleDefs" both consume.

   Deriving them from the frame rather than writing them out a second time is the point: the tensor
   part resolves its momenta through the frame while the scalar dressings go through AngleDefs, and
   if the two disagree the kernel is silently wrong at every momentum -- nothing downstream compares
   them. *)

frameSpatialCosines[frame_, ql_, qs_List, p_, l1_] :=
  MapIndexed[
    Symbol["cosl1p" <> ToString[First[#2]]] ->
      Simplify[Rest[frame[ql]] . Rest[frame[#1]] / (l1 p)]&,
    qs];
