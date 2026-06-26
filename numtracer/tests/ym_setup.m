fields= <|
"Commuting"-> {A[p,{v, c}]},
"Grassmann"->{{cb[p,{c}],c[p,{c}]}}
|>;


truncation=<|
GammaN->{{A,A},{A,A,A},{A,A,A,A},{A,cb,c},{cb,c}},
Propagator->{{A,A},{cb,c}},Rdot->{{A,A},{cb,c}},
S->{{A,A},{A,A,A},{A,A,A,A},{cb,c},{cb,c,A}},
Field->{{}}
|>;


bases=<|
GammaN->{{A,A}->{"AA",1},{A,A,A}->"AAAClass",{A,A,A,A}->"AAAAClass",{A,cb,c}->{"Acbc",1},{cb,c}->"cbc"},
S->{{A,A}->{"AA",1},{A,A,A}->"AAAClass",{A,A,A,A}->"AAAAClass",{A,cb,c}->{"Acbc",1},{cb,c}->"cbc"},
Propagator->{{A,A}->{"AA",1},{cb,c}->"cbc"},
Rdot->{{A,A}->{"AA",1},{cb,c}->"cbc"}
|>;


diagramStyling=<|"Styles"->{A->{Orange},c->{Black,Dashed}}|>;
FSetTexStyles[cb->"\\bar{c}"];


Setup=<|
"FieldSpace"->fields,
"Truncation"->truncation,
"FeynmanRules"->bases,
"DiagramStyling"->diagramStyling
|>;
FSetGlobalSetup[Setup];


SP3Patt[p1e_,p2e_,p3e_]:={Sqrt[(sp[p1,p1]+sp[p2,p2]+sp[p3,p3])/3]}/.{p1:>p1e,p2:>p2e,p3:>p3e}//UseLorentzLinearity//FullSimplify;
SP4Patt[p1e_,p2e_,p3e_,p4e_]:={Sqrt[(sp[p1,p1]+sp[p2,p2]+sp[p3,p3]+sp[p4,p4])/4]}/.{p1:>p1e,p2:>p2e,p3:>p3e,p4:>p4e}//UseLorentzLinearity//FullSimplify;


dressingRules=ReplaceRepeated[#,{
dressing[GammaN,{cb,c},1,{p1_,p2_}]:>-Zc[Sqrt[sp[p2,p2]]]sp[p2,p2],
dressing[GammaN,{A,A},1,{p1_,p2_}]:>ZA[Sqrt[sp[p2,p2]]]sp[p2,p2],

dressing[InverseProp,{cb,c},1,{p1_,p2_}]:>-(Zc[Sqrt[sp[p2,p2]]]sp[p2,p2]+RB[k^2,sp[p2,p2]]Zc[k]),
dressing[InverseProp,{A,A},1,{p1_,p2_}]:>ZA[Sqrt[sp[p2,p2]]]sp[p2,p2]+RB[k^2,sp[p2,p2]]ZA[evP],

dressing[GammaN,{A,cb,c},1,{p1_,p2_,p3_}]:>ZAcbc[p1,p2], 
dressing[GammaN,{A,A,A},1,{p1_,p2_,p3_}]:>ZA3[p1,p2], 
dressing[GammaN,{A,A,A,A},1,{p1_,p2_,p3_,p4_}]:>ZA4[p1,p2,p3] ,

ZAcbc[p1_,p2_]:>ZAcbc@@SP3Patt[p1,p2,-p1-p2],
ZA3[p1_,p2_]:>ZA3@@SP3Patt[p1,p2,-p1-p2],
ZA4[p1_,p2_,p3_]:>ZA4@@SP4Patt[p1,p2,p3,-p1-p2-p3],

nZA->6,
evP:>(k^nZA+1)^(1/nZA),
devP:>k^(-1+nZA) (1+k^nZA)^(-1+1/nZA),
dressing[Rdot,{A,A},1,{p1_,p2_}]:>ZA[evP]RBdot[k^2,sp[p2,p2]]+RB[k^2,sp[p2,p2]](dtZA[evP]+k*devP*(ZA[1.02evP]-ZA[evP])/(0.02*evP)),
dressing[Rdot,{cb,c},1,{p1_,p2_}]:>Zc[k]RBdot[k^2,sp[p2,p2]]+RB[k^2,sp[p2,p2]](dtZc[k]+k (Zc[1.02*k]-Zc[k])/(0.02*k))
}]&;

FSetSymmetricDressing[GammaN,{A,A}]


PropParam[expr_]:=UseLorentzLinearity[expr]//.{
lf1->l1,(*We don't care about this in vacuum*)
sp[p1,p1]->p^2,sp[l1,l1]->l1^2,
sp[l1,p1]->l1 p cos[p,l1],
sp[p1,l1]->l1 p cos[p,l1],
Sqrt[a_^2]:>a,(a_^2)^(n_/2):>a^n,
cos[l1,p]:>cos1
};

SP3FormRule=FMakeSPFormRule[{l1,lf1},p,{p1,p2,p3}];
SP4FormRule=FMakeSPFormRule[{l1,lf1},p,{p1,p2,p3,p4}];
SPParam[expr_]:=UseLorentzLinearity[expr]//.{
lf1->l1,(*We don't care about this in vacuum*)
sp[p,p]->p^2,sp[l1,l1]->l1^2,
sp[l1,p1]->p l1 cos[l1,p1],
sp[l1,p2]->p l1 cos[l1,p2],
sp[l1,p3]->p l1 cos[l1,p3],
sp[l1,p4]->p l1 cos[l1,p4],

Sqrt[a_^2]:>a,(a_^2)^(n_/2):>a^n,(a_^2)^(n_/2):>a^n,Power[Power[l1_,2],Rational[n_,2]]:>l1^n,
cos[l1,p1]:>cosl1p1,
cos[l1,p2]:>cosl1p2,
cos[l1,p3]:>cosl1p3,
cos[l1,p4]:>cosl1p4
};

SetNc[3]
$Assumptions=k>0&&p>0&&l1>0&&-1<cos1<1&&-1<cos2<1&&-1<cos3<1;


(* ::Input::Initialization:: *)
interpolatorType="SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";

(* Codegen opt-in for the TESTS: emit kernels against the standalone shim (DiFfRG namespace +
   complex/compute API) instead of NumTracer's self-contained default, so the generated kernels
   match the copied FormTracer oracles. Applied via SetOptions in each generator AFTER NumTracer loads. *)
(* Dressings: the pure-YM generators (gen_ym_a3/a4/flows) don't pass a dressing list and used to
   rely on the codegen default; that default is now {} (general libraries ship no flow dressings),
   so the gauge set is restored here for the YM tests. *)
$ntShimOpts={"RuntimeInclude"->"shim.hpp","ExtraIncludes"->{},
  "KernelNamespace"->"DiFfRG","SupportNamespace"->"DiFfRG","DressingType"->interpolatorType,
  "Dressings"->{ZA3,ZAcbc,ZA4,dtZc,Zc,dtZA,ZA}};

kernelParameterList={
<|"Name"->"k","Type"->"double"|>,
(*strong couplings*)
<|"Name"->"ZA3","Type"->interpolatorType,"Const"->True,"Reference"->True|>,
<|"Name"->"ZAcbc","Type"->interpolatorType,"Const"->True,"Reference"->True|>,
<|"Name"->"ZA4","Type"->interpolatorType,"Const"->True,"Reference"->True|>,
(*ghost propagator*)
<|"Name"->"dtZc","Type"->interpolatorType,"Const"->True,"Reference"->True|>,
<|"Name"->"Zc","Type"->interpolatorType,"Const"->True,"Reference"->True|>,
(*glue propagator*)
<|"Name"->"dtZA","Type"->interpolatorType,"Const"->True,"Reference"->True|>,
<|"Name"->"ZA","Type"->interpolatorType,"Const"->True,"Reference"->True|>
};

SP4Defs=DeclareSymmetricPoints4DP4[l1,p,{p1,p2,p3,p4}];
SP3Defs=DeclareSymmetricPoints4DP3[l1,p,{p1,p2,p3}];


FSetRegisterSize[64]
