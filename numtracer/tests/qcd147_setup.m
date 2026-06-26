(* ::Package:: *)

(* qcd147_setup.m — TEST fixture, a sibling of qcd_setup.m that uses the FULL quark-gluon
   vertex tensor basis {1,4,7} ("AqbqDirect147") for {A, qb, q} instead of the struct-1-only
   "AqbqDirect1". Everything else (field space, propagators, dressing rules, frames) is identical
   to qcd_setup.m. This is the dense `1/4/7` case the numeric backend must handle; qcd_setup.m is
   kept untouched as the struct-1 regression guard. NOT part of the NumTracer frontend. *)

fields = <|
  "Commuting" -> {A[p, {v, c}]},
  "Grassmann" -> {{cb[p, {c}], c[p, {c}]}, {qb[p, {d, A, F}], q[p, {d, A, F}]}}
|>;

TBRestrictBasis["AqbqDirect", "AqbqDirect147", {1, 4, 7}];
TBRestrictBasis["AqbqDirect", "AqbqDirect1", {1}];

truncation = <|
  GammaN -> {{A, A}, {A, A, A}, {A, A, A, A}, {A, cb, c}, {cb, c}, {A, qb, q}, {qb, q}},
  Propagator -> {{A, A}, {cb, c}, {qb, q}},
  R -> {{A, A}, {cb, c}, {qb, q}},
  Rdot -> {{A, A}, {cb, c}, {qb, q}},
  Field -> {{}}
|>;

bases = <|
  GammaN -> {{A, A} -> {"AA", 1}, {A, A, A} -> "AAAClass", {A, A, A, A} -> "AAAAClass",
     {A, cb, c} -> {"Acbc", 1}, {cb, c} -> "cbc", {A, qb, q} -> "AqbqDirect147", {qb, q} :> "qbq"},
  Propagator -> {{A, A} -> {"AA", 1}, {cb, c} -> "cbc", {qb, q} -> "qbq"},
  R -> {{A, A} -> {"AA", 1}, {cb, c} -> "cbc", {qb, q} -> {"qbq", 1}},
  Rdot -> {{A, A} -> {"AA", 1}, {cb, c} -> "cbc", {qb, q} -> {"qbq", 1}}
|>;

Setup = <|
  "FieldSpace" -> fields,
  "Truncation" -> truncation,
  "FeynmanRules" -> bases,
  "DiagramStyling" -> <|"Styles" -> {A -> {Orange}, c -> {Black, Dashed}, q -> {Black}}|>
|>;
FSetGlobalSetup[Setup];

SP3Patt[p1e_, p2e_, p3e_] := {Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3])/3]} /.
    {p1 :> p1e, p2 :> p2e, p3 :> p3e} // UseLorentzLinearity // FullSimplify;
SP4Patt[p1e_, p2e_, p3e_, p4e_] := {Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3] + sp[p4, p4])/4]} /.
    {p1 :> p1e, p2 :> p2e, p3 :> p3e, p4 :> p4e} // UseLorentzLinearity // FullSimplify;

dressingRules[expr_] := Module[{tmp},
  tmp = ReplaceRepeated[expr, {
    (* Propagators *)
    dressing[GammaN, {cb, c}, 1, {p1_, p2_}] :> -Zc[Sqrt[sp[p2, p2]]] sp[p2, p2],
    dressing[GammaN, {A, A}, 1, {p1_, p2_}] :> ZA[Sqrt[sp[p2, p2]]] sp[p2, p2],
    dressing[GammaN, {qb, q}, 1, {p1_, p2_}] :> -Zq[Sqrt[sp[p2, p2]]],
    dressing[GammaN, {qb, q}, 2, {p1_, p2_}] :> -Mq[Sqrt[sp[p2, p2]]],
    (* Inverse propagators (contain the regulator) *)
    dressing[InverseProp, {cb, c}, 1, {p1_, p2_}] :> -(Zc[Sqrt[sp[p2, p2]]] sp[p2, p2] + RB[k^2, sp[p2, p2]] Zc[k]),
    dressing[InverseProp, {A, A}, 1, {p1_, p2_}] :> ZA[Sqrt[sp[p2, p2]]] sp[p2, p2] + RB[k^2, sp[p2, p2]] ZA[evP],
    dressing[InverseProp, {qb, q}, 1, {p1_, p2_}] :> -(Zq[Sqrt[sp[p2, p2]]] + (Zq[k] RF[k^2, sp[p2, p2]])/Sqrt[sp[p2, p2]]),
    dressing[InverseProp, {qb, q}, 2, {p1_, p2_}] :> -Mq[Sqrt[sp[p2, p2]]],
    (* Strong couplings *)
    dressing[GammaN, {A, cb, c}, 1, {p1_, p2_, p3_}] :> ZAcbc[p1, p2],
    dressing[GammaN, {A, A, A}, 1, {p1_, p2_, p3_}] :> ZA3[p1, p2],
    dressing[GammaN, {A, A, A, A}, 1, {p1_, p2_, p3_, p4_}] :> ZA4[p1, p2, p3],
    dressing[GammaN, {A, qb, q}, 1, {p1_, p2_, p3_}] :> ZAqbq1[p1, p2],
    dressing[GammaN, {A, qb, q}, 2, {p1_, p2_, p3_}] :> ZAqbq4[p1, p2],
    dressing[GammaN, {A, qb, q}, 3, {p1_, p2_, p3_}] :> ZAqbq7[p1, p2],
    (* Parametrizations *)
    ZAcbc[p1_, p2_] :> ZAcbc @@ SP3Patt[p1, p2, -p1 - p2],
    ZA3[p1_, p2_] :> ZA3 @@ SP3Patt[p1, p2, -p1 - p2],
    ZA4[p1_, p2_, p3_] :> ZA4 @@ SP4Patt[p1, p2, p3, -p1 - p2 - p3],
    ZAqbq1[p1_, p2_] :> ZAqbq1 @@ SP3Patt[p1, p2, -p1 - p2],
    ZAqbq4[p1_, p2_] :> ZAqbq4 @@ SP3Patt[p1, p2, -p1 - p2],
    ZAqbq7[p1_, p2_] :> ZAqbq7 @@ SP3Patt[p1, p2, -p1 - p2],
    nZA -> 6,
    evP :> (k^nZA + 1)^(1/nZA),
    devP :> k^(-1 + nZA) (1 + k^nZA)^(-1 + 1/nZA),
    (* Regulator derivatives *)
    dressing[Rdot, {A, A}, 1, {p1_, p2_}] :> ZA[evP] RBdot[k^2, sp[p2, p2]] +
       RB[k^2, sp[p2, p2]] (dtZA[evP] + k*devP*(ZA[1.02 evP] - ZA[evP])/(0.02*evP)),
    dressing[Rdot, {cb, c}, 1, {p1_, p2_}] :> Zc[k] RBdot[k^2, sp[p2, p2]] +
       RB[k^2, sp[p2, p2]] (dtZc[k] + k (Zc[1.02*k] - Zc[k])/(0.02*k)),
    (* Minus sign here due to the vertex basis: it's pdash[p1,d1,d2]! *)
    dressing[Rdot, {qb, q}, 1, {p1_, p2_}] :> -(Zq[k] RFdot[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]] +
       RF[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]] (dtZq[k] + k (Zq[1.02*k] - Zq[k])/(0.02*k))),
    (* Regulators *)
    dressing[R, {qb, q}, 1, {p1_, p2_}] :> -(Zq[k] RF[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]])
  }];
  UseLorentzLinearity[tmp] //. {lf1 -> l1} // Return;
];

FSetSymmetricDressing[GammaN, {A, A}];
FSetSymmetricDressing[GammaN, {cb, c}];
FSetSymmetricDressing[GammaN, {qb, q}];

(* vacuum 1-angle scalar reduction: sp/cos -> {l1, p, cos1} *)
PropParam[expr_] := UseLorentzLinearity[expr] //. {
  lf1 -> l1,
  sp[p1, p1] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> l1 p cos[l1, p],
  sp[p1, l1] -> l1 p cos[l1, p],
  Sqrt[a_^2] :> a, (a_^2)^(n_/2) :> a^n,
  cos[l1, p] :> cos1
};

(* vacuum multi-angle scalar reduction for 3-/4-point vertices: sp/cos -> {l1, p, cosl1p_i}.
   The loop-external cosines cosl1p{1..4} are mapped to the symmetric-point angles (cos1,cos2[,phi])
   by the generator script, matching the frame (sp3Frame/sp4Frame) the tensor part uses. *)
SPParam[expr_] := UseLorentzLinearity[expr] //. {
  lf1 -> l1,
  sp[p, p] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> p l1 cos[l1, p1], sp[l1, p2] -> p l1 cos[l1, p2],
  sp[l1, p3] -> p l1 cos[l1, p3], sp[l1, p4] -> p l1 cos[l1, p4],
  Sqrt[a_^2] :> a, (a_^2)^(n_/2) :> a^n, Power[Power[l1_, 2], Rational[n_, 2]] :> l1^n,
  cos[l1, p1] :> cosl1p1, cos[l1, p2] :> cosl1p2,
  cos[l1, p3] :> cosl1p3, cos[l1, p4] :> cosl1p4
};

(* Compile-time SU(N) sizes: colour SU(3), flavour SU(2) (Nf=2). NumTracer bakes these as
   integer constants (a closed quark loop folds to Nf=2); they are NOT runtime kernel parameters. *)
SetNc[3];
SetNf[2];
$Assumptions = k > 0 && p > 0 && l1 > 0 && -1 < cos1 < 1 && -1 < cos2 < 1;

(* Codegen opt-in for the TESTS: emit kernels against the standalone shim (which supplies the
   `DiFfRG` namespace + complex/compute API) instead of NumTracer's self-contained default, so
   the generated kernels match the copied FormTracer oracles byte-for-byte. Each generator script
   applies it via `SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts]` AFTER NumTracer loads. *)
$ntShimOpts = {"RuntimeInclude" -> "shim.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG",
  "DressingType" -> "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>"};
