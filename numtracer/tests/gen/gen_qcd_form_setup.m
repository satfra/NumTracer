(* Shared setup for the QCD FormTracer (FORM-codegen) clean-cache TIMING scripts. Loads the QCD
   field space / dressingRules / PropParam / SPParam from qcd_setup.m, then (a) restores a SYMBOLIC
   Nf — the FORM oracle keeps Nf a runtime kernel parameter exactly as QCD_vacuum_base/QCD.m does
   (it never calls SetNf), so we undo qcd_setup.m's SetNf[2]; and (b) defines the FORM-codegen
   extras (kernelParameterList / SP*FormRule / SP*Defs / interpolatorType) that qcd_setup.m omits,
   lifted verbatim from QCD.m lines 186-214. These scripts use DiFfRG's MakeKernel/FormTrace, so
   NumTracer.m is intentionally NOT loaded. *)

Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
SetNf[];   (* symbolic Nf: the FORM kernel takes Nf as a runtime double, like QCD.m *)

interpolatorType = "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";

kernelParameterList = {
  <|"Name" -> "k", "Type" -> "double"|>,
  <|"Name" -> "Nf", "Type" -> "double"|>,
  <|"Name" -> "ZA3", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZAcbc", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZA4", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
  <|"Name" -> "ZAqbq1", "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
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
