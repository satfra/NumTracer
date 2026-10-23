#!/usr/bin/env wolframscript
(* Unit gate for the loop-invariant hoisting of finite-T kernels ("Hoist" -> True): ntHoistKernelText
   and its helpers in mathematica/CodegenKernel.m.

   The transform works on the lowered C++ text. Every `const T name = expr;` line that depends neither
   on an excluded parameter nor on a line that does moves into pre(); the rest stays in post(). DiFfRG
   evaluates pre() once per (radial node, frequency node) and post() at every quadrature point, so a
   line wrongly moved into pre() reads a stale or zero value, silently. Pinned here:

     1. parameter lists split at top-level commas only (template arguments carry commas);
     2. the excluded parameters: box variables after the radial one, the external coordinates between
        the frequency and the first runtime parameter, and "PointArguments";
     3. dependence is transitive, and a name declared by any other line (the fenv array) counts as
        post-only;
     4. pre_reads flags exactly the parameters pre() reads, in kernel() order;
        pre() returns only the values post() reads;
     5. every definition of a complex flow's #if ladder is transformed;
     6. a kernel without the frequency parameter is refused (None), not half-transformed;
     7. the real emission shapes: shared lookups carry [[maybe_unused]], and k, which reaches
        ntHoistedKernelFn outside ScalarParams, is a runtime parameter, not an external coordinate;
     8. a vacuum flow (no frequency) hoists over the radial variable alone.

   Run standalone:  wolfram -script numtracer/tests/test_hoist_split.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

$nExpected = 36;
$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);

runtime = {"k", "T", "d1V", "ZQ", "ZPi"};

Print["=== 1. top-level parameter split ==="];
check["template commas stay inside a parameter",
  ntSplitTopLevel["const double& l1, const Foo<double, Bar<int, 2>>& ZQ, const auto& d1V"] ===
    {"const double& l1", "const Foo<double, Bar<int, 2>>& ZQ", "const auto& d1V"}];

Print["=== 2. excluded parameters ==="];
check["map kernel: angle and external coordinates",
  ntHoistExcluded[{"l1", "cos1", "f0", "p0", "p", "k", "T", "d1V", "ZQ"}, "f0", runtime, {}] === {"cos1", "p0", "p"}];
check["potential kernel: only the point arguments",
  ntHoistExcluded[{"l1", "f0", "k", "T", "d1V", "ZQ"}, "f0", runtime, {"d1V"}] === {"d1V"}];
check["no frequency parameter: None",
  ntHoistExcluded[{"l1", "cos1", "k"}, "f0", runtime, {}] === None];

kernel = "static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& f0, const double& p, const double& k, const double& T, const auto& d1V, const auto& ZQ, const auto& ZPi)
    {
      using namespace DiFfRG;
      double fenv[3];
      const double dr_0 = ZQ(f0 + 3.14 * T, l1);
      const double dr_1 = ZQ(f0, sqrt(l1 * l1 + p * p + 2. * l1 * p * cos1));
      DiFfRG::ns::fill(fenv, l1, cos1, dr_0, dr_1);
      const auto _a = powr<-1>(dr_0 + k);
      const auto _b = dr_1 * _a;
      const auto _c = DiFfRG::ns::tr0(fenv) * _a;
      const auto _d = _a + d1V;
      return _b + _c + _d;
    }";

Print["=== 3. a map kernel ==="];
res = ntHoistKernelText[kernel, "f0", runtime, {"d1V"}];
pre = StringCases[res, Shortest["auto pre(" ~~ __ ~~ "\n}"]];
post = StringCases[res, Shortest["auto post(" ~~ __ ~~ "\n}"]];
check["one pre() and one post()", Length[pre] === 1 && Length[post] === 1];
check["loop-only lines go to pre()", StringContainsQ[First[pre], "const double dr_0"] && StringContainsQ[First[pre], "const auto _a"]];
check["a line reading the angle stays in post()", !StringContainsQ[First[pre], "dr_1"] && StringContainsQ[First[post], "const double dr_1"]];
check["dependence is transitive", !StringContainsQ[First[pre], "_b ="] && StringContainsQ[First[post], "const auto _b"]];
check["a line reading the fenv array stays in post()", !StringContainsQ[First[pre], "_c ="] && StringContainsQ[First[post], "const auto _c"]];
check["a line reading a point argument stays in post()", !StringContainsQ[First[pre], "_d ="] && StringContainsQ[First[post], "const auto _d"]];
(* both are read by post(): dr_0 by the fill() call *)
check["pre() returns what post() reads", StringContainsQ[First[pre], "return DiFfRG::device::make_tuple(dr_0, _a);"]];
check["post() binds them in order",
  StringContainsQ[First[post], "const auto &dr_0 = DiFfRG::device::get<0>(P_);"] &&
  StringContainsQ[First[post], "const auto &_a = DiFfRG::device::get<1>(P_);"]];
check["the return stays in post()", StringContainsQ[First[post], "return _b + _c + _d;"]];
(* l1 f0 T k ZQ read; cos1 p d1V ZPi not *)
check["pre_reads in kernel() order",
  StringContainsQ[res, "static constexpr std::array<bool, 9> pre_reads{true, false, true, false, true, true, false, true, false};"]];
check["the original kernel() is kept", StringContainsQ[res, kernel]];

Print["=== 4. a complex flow's #if ladder ==="];
ladder = "#if X == 2\n" <> kernel <> "\n#else\n" <> kernel <> "\n#endif";
res2 = ntHoistKernelText[ladder, "f0", runtime, {"d1V"}];
check["both definitions transformed", StringCount[res2, "auto pre("] === 2 && StringCount[res2, "auto post("] === 2];
check["the preprocessor lines survive", StringCount[res2, "#if X == 2"] === 1 && StringCount[res2, "#endif"] === 1];

Print["=== 5. refusal ==="];
check["no frequency parameter: None", ntHoistKernelText[StringReplace[kernel, "f0" -> "q0"], "f0", runtime, {}] === None];

Print["=== 6. what pre() returns ==="];
inner = "static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& f0, const double& p, const double& k, const double& T, const auto& ZQ)
    {
      const double dr_0 = f0 + 3.14 * T;
      const auto _a = ZQ(dr_0, l1);
      const auto _u = l1 * T;
      return _a * cos1 * p;
    }";
res4 = ntHoistKernelText[inner, "f0", {"k", "T", "ZQ"}, {}];
pre4 = First @ StringCases[res4, Shortest["auto pre(" ~~ __ ~~ "\n}"]];
check["a value read only inside pre() is not returned", StringContainsQ[pre4, "return DiFfRG::device::make_tuple(_a);"]];
check["a value nothing reads is dropped", !StringContainsQ[pre4, "_u ="]];
check["post() reads nothing from pre(): not hoisted",
  ntHoistKernelText[StringReplace[inner, "return _a * cos1 * p;" -> "return cos1 * p;"], "f0", {"k", "T", "ZQ"}, {}] === None];

Print["=== 7. emission shapes ==="];
emitted = "static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& f0, const double& p, const double& k, const double& T, const auto& d1V, const auto& ZQ)
    {
      using namespace DiFfRG;
      [[maybe_unused]] const auto ntgix0 = ZQ.index(f0, l1);
      [[maybe_unused]] const auto ntglk0 = ZQ.at(ntgix0);
      const auto _r = RB(powr<2>(k), powr<2>(l1));
      return _r * ntglk0 * (d1V + cos1 * p);
    }";
res3 = First @ ntHoistedKernelFn[emitted, <|"Hoist" -> True, "ScalarParams" -> {"T", "d1V"}, "Dressings" -> {"ZQ"}, "PointArguments" -> {}|>, False, f0];
pre3 = First @ StringCases[res3, Shortest["auto pre(" ~~ __ ~~ "\n}"]];
check["[[maybe_unused]] lookups go to pre()", StringContainsQ[pre3, "ntgix0 = ZQ.index"] && StringContainsQ[pre3, "ntglk0 = ZQ.at"]];
check["a line reading k goes to pre()", StringContainsQ[pre3, "const auto _r = RB"]];
check["pre_reads flags k", StringContainsQ[res3, "pre_reads{true, false, true, false, true, false, false, true}"]];

Print["=== 8. vacuum flows ==="];
check["vacuum map kernel: angles and external coordinate",
  ntHoistExcluded[{"l1", "cos1", "cos2", "p", "k", "T", "d1V", "ZQ"}, None, runtime, {}] === {"cos1", "cos2", "p"}];
check["vacuum potential kernel: only the point arguments",
  ntHoistExcluded[{"l1", "k", "T", "d1V", "ZQ"}, None, runtime, {"d1V"}] === {"d1V"}];
check["vacuum: no radial parameter: None", ntHoistExcluded[{"k", "T", "d1V"}, None, runtime, {}] === None];
vac = "static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const double& T, const auto& d1V, const auto& ZQ)
    {
      using namespace DiFfRG;
      const double dr_0 = ZQ(l1);
      const double dr_1 = ZQ(sqrt(l1 * l1 + p * p + 2. * l1 * p * cos1));
      const auto _a = powr<-1>(dr_0 + k);
      const auto _b = dr_1 * _a + d1V;
      return _b;
    }";
{resV, traitV} = ntHoistedKernelFn[vac, <|"Hoist" -> True, "ScalarParams" -> {"T", "d1V"}, "Dressings" -> {"ZQ"}, "PointArguments" -> {"d1V"}|>, False, None];
preV = First @ StringCases[resV, Shortest["auto pre(" ~~ __ ~~ "\n}"]];
postV = First @ StringCases[resV, Shortest["auto post(" ~~ __ ~~ "\n}"]];
check["vacuum: radial-only lines go to pre()", StringContainsQ[preV, "const double dr_0"] && StringContainsQ[preV, "const auto _a"]];
check["vacuum: a line reading an angle or p stays in post()", !StringContainsQ[preV, "dr_1"] && StringContainsQ[postV, "const double dr_1"]];
check["vacuum: a line reading a point argument stays in post()", StringContainsQ[postV, "const auto _b"]];
(* l1 k ZQ read; cos1 p T d1V not *)
check["vacuum: pre_reads in kernel() order", StringContainsQ[resV, "pre_reads{true, false, false, true, false, false, true}"]];
check["vacuum: the has_pre trait", traitV === {"static constexpr bool has_pre = true;"}];

Print["=== explicit template parameters ==="];
(* the form ntExplicitTemplates gives the kernel: a `template<...>` line above the definition, here
   inside an #if branch. pre() must keep the header, post() must add P_'s type to it, and no
   abbreviated `auto` parameter may appear (nvcc 12.x cudafe++ crashes on many of them). *)
tmpl = "#if 1
    template<typename T7, typename T8>
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const double& T, const double& d1V, const T7& ZQ, const T8& ZPi)
    {
      const double dr_0 = ZQ(l1);
      const auto _a = powr<-1>(dr_0 + k);
      const auto _b = cos1 * _a;
      return _b;
    }
#endif";
{resT, traitT} = ntHoistedKernelFn[tmpl, <|"Hoist" -> True, "ScalarParams" -> {"T", "d1V"}, "Dressings" -> {"ZQ", "ZPi"}, "PointArguments" -> {}|>, False, None];
check["template: pre() keeps the header",
  StringContainsQ[resT, "template<typename T7, typename T8>\n    static KOKKOS_INLINE_FUNCTION auto pre("]];
check["template: post() adds P_'s type to the header",
  StringContainsQ[resT, "template<typename TP_, typename T7, typename T8>\n    static KOKKOS_INLINE_FUNCTION auto post(const TP_ &P_, "]];
check["template: no abbreviated auto parameter", StringFreeQ[resT, "(const auto" | ", const auto"]];
check["no template header: post() gets its own",
  StringContainsQ[resV, "template<typename TP_>\n"] && StringContainsQ[resV, "auto post(const TP_ &P_, "]];

If[nRan =!= $nExpected,
  Print["  FAIL  check count is ", nRan, ", expected ", $nExpected,
        " — a check was added/removed (update $nExpected) or an Abort escaped mid-file"];
  $fail++,
  Print["  PASS  all ", $nExpected, " checks ran (no early Abort)"]];
Print[""];
If[TrueQ[$fail > 0],
  Print["FAILED: ", $fail, " of ", nRan, " check(s)"]; End[]; Exit[1],
  Print["all ", nRan, " checks passed"]; End[]; Exit[0]];
