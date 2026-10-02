(* ::Package:: *)

(* qcd_lambda3d_setup.m — TEST fixture: the four-quark flow at TWO-external-momentum kinematics.
   Cloned from qcd_aaqbq_setup.m and trimmed, for the same reason that file was cloned from
   qcd_setup.m: the frozen setups guard committed kernels, and changing a truncation changes the
   Wetterich derivative of every flow generated from it.

   WHY THIS FIXTURE EXISTS. It is the only thing in the suite that exercises the GENERAL-FRAME
   (polyFrameSpec) code path: every other committed flow satisfies unitLoopOkQ (externals depend on
   the single magnitude `p`), so the unit-loop parametrisation with its unit-constraint groups is
   taken and `reduce_units` collapses sin^2 -> 1-cos^2 and the bare-loop denominator to l1^2. A flow
   whose externals are parametrised by shape coordinates (here S0,S1,SPhi — the with_mesons
   lambda1L3D class) falls off that path, and before this fixture NOTHING in the suite covered what
   happens there: minted radical symbols (ntSin$/ntRad$), a wider Poly symbol space, momentum bases
   appearing in exact +/- pairs, and dressed+complex traces of real size. Those are exactly the
   mechanisms the frame/sign levers work on, so they need a guard.

   CONTENT. The derivative is the genuine four-quark object {qb,q,qb,q} (the lambda1 channel),
   projected on channel 1 of the Fierz-complete Nf=2 basis, at the two-momentum configuration
   (p2,-p2,-p3,p3) of with_mesons' lambda1tch3D — gluon-exchange boxes through the AqbqDirect1
   vertex. No mesons: the meson sector only adds more dressing atoms, while every mechanism under
   test (frame fallback, +/- momentum routing, dress atoms, complex traces, sub-term dedup) is
   already present with the gauge sector alone, and this keeps the fixture cheap. *)

(* Pin FunKit's backend BEFORE any derivative/trace is taken: FRoute is leg-order sensitive, so
   the backend choice silently changes the loop-momentum routing and thus every frozen test
   integrand. See backend_pin.m. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "backend_pin.m"}]];

fields = <|
  "Commuting" -> {A[p, {v, c}]},
  "Grassmann" -> {{qb[p, {d, A, F}], q[p, {d, A, F}]}}
|>;

TBRestrictBasis["AqbqDirect", "AqbqDirect1", {1}];

truncation = <|
  GammaN -> {{A, A}, {A, qb, q}, {qb, q}},
  Propagator -> {{A, A}, {qb, q}},
  R -> {{A, A}, {qb, q}},
  Rdot -> {{A, A}, {qb, q}},
  Field -> {{}}
|>;

bases = <|
  GammaN -> {{A, A} -> {"AA", 1}, {A, qb, q} -> "AqbqDirect1", {qb, q} :> "qbq"},
  Propagator -> {{A, A} -> {"AA", 1}, {qb, q} -> "qbq"},
  R -> {{A, A} -> {"AA", 1}, {qb, q} -> {"qbq", 1}},
  Rdot -> {{A, A} -> {"AA", 1}, {qb, q} -> {"qbq", 1}}
|>;

Setup = <|
  "FieldSpace" -> fields,
  "Truncation" -> truncation,
  "FeynmanRules" -> bases,
  "DiagramStyling" -> <|"Styles" -> {A -> {Orange}, q -> {Black}}|>
|>;
FSetGlobalSetup[Setup];

SP3Patt[p1e_, p2e_, p3e_] := {Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3])/3]} /.
    {p1 :> p1e, p2 :> p2e, p3 :> p3e} // UseLorentzLinearity // FullSimplify;

dressingRules[expr_] := Module[{tmp},
  tmp = ReplaceRepeated[expr, {
    (* Propagators *)
    dressing[GammaN, {A, A}, 1, {p1_, p2_}] :> ZA[Sqrt[sp[p2, p2]]] sp[p2, p2],
    dressing[GammaN, {qb, q}, 1, {p1_, p2_}] :> -Zq[Sqrt[sp[p2, p2]]],
    dressing[GammaN, {qb, q}, 2, {p1_, p2_}] :> -Mq[Sqrt[sp[p2, p2]]],
    (* Inverse propagators (contain the regulator) *)
    dressing[InverseProp, {A, A}, 1, {p1_, p2_}] :> ZA[Sqrt[sp[p2, p2]]] sp[p2, p2] + RB[k^2, sp[p2, p2]] ZA[evP],
    dressing[InverseProp, {qb, q}, 1, {p1_, p2_}] :> -(Zq[Sqrt[sp[p2, p2]]] + (Zq[k] RF[k^2, sp[p2, p2]])/Sqrt[sp[p2, p2]]),
    dressing[InverseProp, {qb, q}, 2, {p1_, p2_}] :> -Mq[Sqrt[sp[p2, p2]]],
    (* The quark-gluon coupling *)
    dressing[GammaN, {A, qb, q}, 1, {p1_, p2_, p3_}] :> ZAqbq1[p1, p2],
    ZAqbq1[p1_, p2_] :> ZAqbq1 @@ SP3Patt[p1, p2, -p1 - p2],
    nZA -> 6,
    evP :> (k^nZA + 1)^(1/nZA),
    devP :> k^(-1 + nZA) (1 + k^nZA)^(-1 + 1/nZA),
    (* Regulator derivatives *)
    dressing[Rdot, {A, A}, 1, {p1_, p2_}] :> ZA[evP] RBdot[k^2, sp[p2, p2]] +
       RB[k^2, sp[p2, p2]] (dtZA[evP] + k*devP*(ZA[1.02 evP] - ZA[evP])/(0.02*evP)),
    (* Minus sign here due to the vertex basis: it's pdash[p1,d1,d2]! *)
    dressing[Rdot, {qb, q}, 1, {p1_, p2_}] :> -(Zq[k] RFdot[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]] +
       RF[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]] (dtZq[k] + k (Zq[1.02*k] - Zq[k])/(0.02*k))),
    (* Regulators *)
    dressing[R, {qb, q}, 1, {p1_, p2_}] :> -(Zq[k] RF[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]])
  }];
  UseLorentzLinearity[tmp] //. {lf1 -> l1} // Return;
];

FSetSymmetricDressing[GammaN, {A, A}];
FSetSymmetricDressing[GammaN, {qb, q}];

(* ===== the (S0,S1,SPhi) frame — copied VERBATIM from with_mesons/QCD.wl =====
   The frame fixes the tensor network; the sp-rule list below reduces the SCALAR part (propagator
   denominators and dressing arguments). The two MUST describe the same kinematics, or every trace
   evaluates its Dirac/Lorentz structure and its dressings at different momenta — silently. The
   frame is built from the very Q/k decomposition TB3PToS0S1SPhi uses, making the two consistent by
   construction (see the QCD.wl comment block for the reflection bug this avoids). *)
frameS0S1Sphi = Module[{aa, ss, Qm, km, zz, Qv, kv, lv},
  aa = S1*Cos[SPhi]; ss = S1*Sin[SPhi];
  Qm = Sqrt[S0^2*(1 - ss)];
  km = Sqrt[(3/4)*S0^2*(1 + ss)];
  zz = aa/Sqrt[1 - ss^2];
  Qv = Qm*{1, 0, 0, 0};
  kv = km*{zz, Sqrt[1 - zz^2], 0, 0};
  lv = l1*{cos1, Sqrt[1 - cos1^2]*cos2, Sqrt[(1 - cos1^2)*(1 - cos2^2)], 0};
  <|p1 -> Qv, p2 -> -kv - Qv/2, p3 -> kv - Qv/2, l1 -> lv|>];
argsS0S1SPhi = {l1, cos1, cos2, S0, S1, SPhi, k};

(* Scalar reduction: value each sp[.,.] through the same TB3PToS0S1SPhi map (sp-only rules, so
   tensor slots keep p1/p2/p3/l1 and stay resolvable by the frame — see QCD.wl for why the map must
   not be applied to a whole trace). sp[l1,l1] is supplied by hand. *)
transf3PTo = TB3PToS0S1SPhi[p1, p2, p3, l1, S0, S1, SPhi];
S0SpRules = Join[{sp[l1, l1] -> l1^2},
  Flatten@Table[With[{v = FullSimplify[transf3PTo[sp[a, b]]]}, {sp[a, b] -> v, sp[b, a] -> v}],
    {a, {p1, p2, p3}}, {b, {p1, p2, p3, l1}}]];
S0Param[expr_] := UseLorentzLinearity[expr] //. Join[{lf1 -> l1}, S0SpRules] //. {
   Sqrt[(a_)^2] :> a, ((a_)^2)^((n_)/2) :> a^n};

(* Compile-time SU(N) sizes: colour SU(3), flavour SU(2) — the Fierz-complete basis needs both. *)
SetNc[3];
SetNf[2];
$Assumptions = k > 0 && l1 > 0 && -1 < cos1 < 1 && -1 < cos2 < 1 &&
  S0 > 0 && 0 < S1 < 1 && -Pi < SPhi < Pi;

$ntShimOpts = {"RuntimeInclude" -> "shim.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG",
  "DressingType" -> "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>"};
