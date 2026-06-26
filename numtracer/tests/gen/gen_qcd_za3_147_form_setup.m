(* FormTracer (FORM-codegen) setup for the QCD three-gluon vertex with the FULL quark-gluon vertex
   {1,4,7} tensor basis — the oracle for the dense `1/4/7` numeric kernel. Sibling of
   gen_qcd_form_setup.m, but loads qcd147_setup.m (so {A,qb,q} uses AqbqDirect147) and carries the
   ZAqbq1/ZAqbq4/ZAqbq7 dressings. Nf is restored symbolic (runtime kernel parameter, like QCD.m). *)

Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
SetNf[];   (* symbolic Nf: the FORM kernel takes Nf as a runtime double, like QCD.m *)

interpolatorType = "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";

kernelParameterList = {
  <|"Name" -> "k", "Type" -> "double"|>,
  <|"Name" -> "Nf", "Type" -> "double"|>,
  <|"Name" -> "ZA3", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZAcbc", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZA4", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZAqbq1", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZAqbq4", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZAqbq7", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "dtZc", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "Zc", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "dtZA", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZA", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "dtZq", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "Zq", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "Mq", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>};

SP3FormRule = FMakeSPFormRule[{l1}, p, {p1, p2, p3}];
SP4FormRule = FMakeSPFormRule[{l1}, p, {p1, p2, p3, p4}];
SP4Defs = DeclareSymmetricPoints4DP4[l1, p, {p1, p2, p3, p4}];
SP3Defs = DeclareSymmetricPoints4DP3[l1, p, {p1, p2, p3}];
FSetRegisterSize[64];
