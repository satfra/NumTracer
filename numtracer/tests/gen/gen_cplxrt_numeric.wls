#!/usr/bin/env wolframscript
(* Codegen gate for FINITE-DENSITY coefficients — a Complex nested below a Power — and for the two
   options that exist because of them: "ComplexRuntimeProjection" and "RealOutput".

   WHAT IS BEING PINNED, and why nothing else pins it.

   Every other complex fixture carries its `I` as a bare multiplicative factor. The re/im split can
   read that off directly, and does. A finite-DENSITY quark propagator cannot be written that way: the
   chemical potential enters as `l0 + I muq` inside a DENOMINATOR, i.e. a Complex under Power[·, -n].
   Substituting I -> ii and reading Coefficient[·, ii, 0|1] is exact only on a POLYNOMIAL in ii; on a
   rational function Mathematica silently power-series expands and returns the leading term, so
       "tr0"/(a + I b)   projects to   ntRe["tr0"]/a
   — muq is deleted from every propagator and the kernel is a wrong number that compiles. ntImagSplitSafeQ
   refuses that shape; "ComplexRuntimeProjection" instead keeps the coefficient intact and emits
   ntRe/ntIm around it, deferring that one projection to generated-code runtime.
   Until this fixture, NO generator anywhere put a Complex under a Power — the option was covered only
   by the symbolic unit test (tests/test_complex_projection.wls §5). This is the end-to-end half:
   it proves the deferred form survives FunKit lowering, the probe, and the C++ compiler, and that it
   still computes the right number.

   THE ORACLE is an identity, not a tolerance: the untouched complex body performs every complex
   multiplication in C++ and is correct by construction, so
       kernel<RePart>  ==  real(kernel<complex>)
   must hold for any integrand. compare_cplxrt_num.cpp forces the two verdicts into two TUs and
   compares them. Same mechanism as the cplxdisc gate.

   EXPECT VERDICT 0. With a genuinely complex denominator Im(integrand) does not vanish, so the probe
   correctly says "complex". That is the point rather than a problem: the grader forces the verdicts it
   needs, and verdict 0 is exactly the case "RealOutput" has to handle. So assert that the bodies were
   BUILT (>= 2 VERDICT branches), never that the committed verdict is 1.

   THE SECOND KERNEL is the "RealOutput" gate. The same ntk is emitted again with that option on,
   which drops the untouched complex body — the expensive one to lower, and the one a real-valued
   consumer could never instantiate. Under verdict 0 it must fall through to the RePart body behind a
   #warning, and must still equal the same complex oracle. That is the whole contract.

   Run from anywhere:  wolfram -script tests/gen/gen_cplxrt_numeric.wls                            *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

(* The finite-T propagator frame: slot 0 is Euclidean time, so the loop carries a genuine l0. Do NOT
   hand-roll this — Codegen*.m picks the component parametrisation via unitLoopMixedOkQ, which needs
   the loop's spatial part to be a unit 3-vector times l1. *)
frame = NumTracer`propFrameFT[Global`p0, Global`p, Global`l0, Global`l1, Global`cos1,
                              Global`pm, Global`lm];

(* muq and Ep are frame-free: they never appear in the frame, only in the coefficient below. Unused
   args are harmless — the emitted fill() marks every one [[maybe_unused]]. The probe draws every arg
   whose name is neither `cos*` nor `phi*` from U(0.25, 3.0), so l0 + I muq + Ep stays away from its
   pole and the sampled kinematics are non-degenerate. *)
args = {Global`l0, Global`l1, Global`cos1, Global`p0, Global`p, Global`muq, Global`Ep};

(* A closed Dirac trace over the two frame momenta, times a fundamental colour loop. tr(T^a T^b T^c)
   f^{abc} is PURELY IMAGINARY, which is what makes the TRACE VALUES themselves complex — a real
   coefficient alone would leave ntIm(tr) identically zero and collapse RePart onto Pure. *)
dirac = ntGamma[m1, s1, s2] ntVec[Global`lm, m1] ntGamma[m2, s2, s1] ntVec[Global`pm, m2];
colT  = ntSUNT[3, a1, c1, c2] ntSUNT[3, a2, c2, c3] ntSUNT[3, a3, c3, c1] ntSUNf[3, a1, a2, a3];

(* THE POINT OF THE FIXTURE. A NEGATIVE INTEGER power of a Complex-bearing base is precisely the shape
   ntImagSplitSafeQ rejects (Times/Plus/non-negative-integer Power are the admissible ones), so this
   is what turns "ComplexRuntimeProjection" on. Squared rather than linear because the option's own
   documentation calls out a denominator raised to a power as the case needing the type-generic powr.
   A scalar factor like this lands untouched in the diagram's "Coeff" via scalarQ. *)
denom = 1/(Global`l0 + I*Global`muq + Global`Ep)^2;

net = denom*dirac*colT + I*dirac*colT;

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];

(* A real coefficient everywhere would silently defeat the fixture: complexQ would go False, the
   bodies would collapse to one, and the gate would pass while testing nothing. *)
If[FreeQ[ntk[[1]]["Diagrams"], Complex],
  Print["FATAL: the fixture is not complex — MakeNTKernel would emit a single body and this gate ",
        "would test nothing. Check that the I and the fundamental-generator colour loop survived."];
  Exit[1]];

(* ... and a coefficient whose Complex is merely a bare factor would exercise the ORDINARY path, not
   the deferral. Assert the nested shape reached the diagrams, or this fixture is just another
   cplxdisc. *)
If[! MemberQ[Cases[ntk[[1]]["Diagrams"], Power[b_, _Integer?Negative] :> b, Infinity],
             _?(! FreeQ[#, Complex] &)],
  Print["FATAL: no Complex under a negative Power survived into the diagrams — the coefficient was ",
        "simplified away, so ComplexRuntimeProjection would never engage and this gate would be a ",
        "duplicate of cplxdisc_num."];
  Exit[1]];

Print["[gen] Numeric backend (ComplexRuntimeProjection) ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_cplxrt_num.cpp"}],
  FileNameJoin[{$gen, "Cplxrt_num_kernel.hh"}],
  FileNameJoin[{$gen, "Cplxrt_num_kernels.hh"}],
  "Name" -> "Cplxrt_num_kernel", "Namespace" -> "cplxrt_num",
  "Dressings" -> {}, "ComplexRuntimeProjection" -> True];

(* The SAME ntk, emitted again with "RealOutput" -> True: two real bodies, no complex body. Separate
   namespace so it gets its own verdict macro and its own probe, and so both can be linked into one
   test binary. *)
Print["[gen] Numeric backend (ComplexRuntimeProjection + RealOutput) ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_cplxrtro_num.cpp"}],
  FileNameJoin[{$gen, "Cplxrtro_num_kernel.hh"}],
  FileNameJoin[{$gen, "Cplxrtro_num_kernels.hh"}],
  "Name" -> "Cplxrtro_num_kernel", "Namespace" -> "cplxrtro_num",
  "Dressings" -> {}, "ComplexRuntimeProjection" -> True, "RealOutput" -> True];

(* The SAME ntk a third time, with "ComplexEndProjection": no symbolic projection at all, one body
   ntRe[integrand] evaluated in complex arithmetic, no probe and no verdict. Its value must equal the
   same complex oracle — compare_cplxrt_num_endproj.cpp. *)
Print["[gen] Numeric backend (ComplexEndProjection + RealOutput) ..."];
Quiet[DeleteFile[FileNameJoin[{$gen, "probe_cplxrtend_num.cpp"}]]];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_cplxrtend_num.cpp"}],
  FileNameJoin[{$gen, "Cplxrtend_num_kernel.hh"}],
  FileNameJoin[{$gen, "Cplxrtend_num_kernels.hh"}],
  "Name" -> "Cplxrtend_num_kernel", "Namespace" -> "cplxrtend_num",
  "Dressings" -> {}, "ComplexRuntimeProjection" -> True, "ComplexEndProjection" -> True,
  "RealOutput" -> True];

Module[{h = Import[FileNameJoin[{$gen, "Cplxrt_num_kernel.hh"}], "Text"],
        hro = Import[FileNameJoin[{$gen, "Cplxrtro_num_kernel.hh"}], "Text"],
        deferred},
(* Both bodies must exist. NOT "the committed verdict is 1" — see the header comment. *)
  If[StringCount[h, "VERDICT"] < 2,
    Print["FATAL: fewer than two VERDICT branches in the emitted kernel — the flow was not treated ",
          "as complex, so neither real projection was built."]; Exit[1]];
  If[StringContainsQ[h, "tr$"],
    Print["FATAL: leaked Unique[\"tr$\"] placeholder(s) in the emitted kernel."]; Exit[1]];

(* THE deferral assertion. ntRe/ntIm applied to a trace call is the ordinary RePart split and proves
   nothing here; what this fixture is for is ntRe/ntIm applied to something that is NOT a trace — the
   surviving complex denominator. Without this the gate would pass on a kernel where the denominator
   had been silently projected (i.e. muq deleted) after all. *)
(* Line-based, NOT a regex on the argument: the deferred expression contains nested parentheses
   (`ntRe(powr<-2>(Ep + l0 + complex<double>(0.,1.) * muq))`), so a `[^)]*` capture stops inside
   complex<double>(0.,1.) and never sees muq. One emitted statement per line makes this exact. *)
  deferred = Select[StringSplit[h, "\n"],
    StringContainsQ[#, "ntRe(" | "ntIm("] && StringContainsQ[#, "muq"] &];
  If[deferred === {},
    Print["FATAL: no ntRe(...)/ntIm(...) in the emitted kernel wraps an expression containing muq. ",
          "The complex denominator was projected symbolically instead of being deferred to runtime, ",
          "which is the exact silent-wrong-kernel bug this fixture exists to catch."];
    Exit[1]];
(* The deferred coefficient must stay FACTORED and complex-typed: an imaginary literal inside a
   negative power. If COEN had expanded it, or the literal had degraded to a real, the runtime split
   would be projecting the wrong thing. This also pins the type-generic powr (upstream issue 5) —
   powr<-2> is being applied to a std::complex here, which the old double-only signature rejected. *)
  If[! AnyTrue[deferred, StringContainsQ[#, "powr<-2>"] && StringContainsQ[#, "complex<double>(0.,1.)"] &],
    Print["FATAL: the deferred coefficient is not `powr<-2>(... complex<double>(0.,1.) ...)`. Either ",
          "the denominator was expanded (defeating the CSE the deferral exists to preserve) or the ",
          "imaginary unit did not survive lowering. Emitted: ", First[deferred]];
    Exit[1]];

(* "RealOutput": the untouched complex body must be GONE. It is the only branch that returns a
   complex type, so its absence is what the emitted class shows. *)
(* Count BODIES, not occurrences of "VERDICT": both spellings mention the macro exactly twice — the
   default as `== 2` / `== 1`, the RealOutput one as `== 2` / the nested `== 0` warning guard — so a
   VERDICT count cannot tell them apart. Three bodies vs two is the actual difference. *)
  If[StringCount[hro, "auto kernel("] =!= 2,
    Print["FATAL: the RealOutput kernel emitted ", StringCount[hro, "auto kernel("],
          " bodies, expected exactly 2 (Pure and RePart, no untouched complex body)."]; Exit[1]];
  If[StringCount[h, "auto kernel("] =!= 3,
    Print["FATAL: the default kernel emitted ", StringCount[h, "auto kernel("],
          " bodies, expected 3 — RealOutput would then not be removing anything."]; Exit[1]];
  If[! StringContainsQ[hro, "warning"],
    Print["FATAL: the RealOutput kernel carries no #warning for the verdict-0 truncation."];
    Exit[1]]];

(* "ComplexEndProjection": ONE body returning ntRe(...) of the whole complex integrand, with no probe and
   no verdict, and the complex denominator left intact (not projected symbolically). *)
Module[{he = Import[FileNameJoin[{$gen, "Cplxrtend_num_kernel.hh"}], "Text"]},
  If[StringCount[he, "auto kernel("] =!= 1,
    Print["FATAL: the ComplexEndProjection kernel emitted ", StringCount[he, "auto kernel("],
          " bodies, expected exactly 1."]; Exit[1]];
  If[StringContainsQ[he, "VERDICT" | "numtrace_verdict.hh"],
    Print["FATAL: the ComplexEndProjection kernel references the probe verdict; it must not need one."];
    Exit[1]];
  If[FileExistsQ[FileNameJoin[{$gen, "probe_cplxrtend_num.cpp"}]],
    Print["FATAL: a probe source was written for the ComplexEndProjection kernel; the mode skips the probe."];
    Exit[1]];
  If[! StringContainsQ[he, "return ntRe("],
    Print["FATAL: the ComplexEndProjection body does not return ntRe(...) of the assembled integrand."];
    Exit[1]];
  If[! AnyTrue[StringSplit[he, "\n"],
        StringContainsQ[#, "powr<-2>"] && StringContainsQ[#, "complex<double>(0.,1.)"] &&
          StringContainsQ[#, "muq"] &],
    Print["FATAL: the complex denominator powr<-2>(... complex<double>(0.,1.) * muq ...) did not survive ",
          "into the ComplexEndProjection kernel; it was projected symbolically."]; Exit[1]]];

(* ... and the mode is refused without "RealOutput": its single body returns a real value. *)
If[! TrueQ[CheckAbort[
      Quiet[MakeNTKernel[ntk,
        FileNameJoin[{$TemporaryDirectory, "gen_cplxrtendbad_num.cpp"}],
        FileNameJoin[{$TemporaryDirectory, "Cplxrtendbad_num_kernel.hh"}],
        FileNameJoin[{$TemporaryDirectory, "Cplxrtendbad_num_kernels.hh"}],
        "Name" -> "Cplxrtendbad_num_kernel", "Namespace" -> "cplxrtendbad_num", "Dressings" -> {},
        "ComplexRuntimeProjection" -> True, "ComplexEndProjection" -> True]]; False,
      True]],
  Print["FATAL: ComplexEndProjection without RealOutput was not refused."]; Exit[1]];

Print["Finite-density (nested-Complex) codegen kernels generated."];
