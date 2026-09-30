#!/usr/bin/env wolframscript
(* Codegen gate for the COMPLEX x MULTI-TRACE intersection — the combination the suite was blind to.

   Two properties have to hold AT ONCE for the real/imaginary projections to be exercised on a
   non-linear integrand, and until this fixture existed no flow had both:

     * MULTI-TRACE: the diagram is DISCONNECTED, so the assembly emits `Π traceRef[...]` and the
       integrand is a POLYNOMIAL, not a linear form, in the trace tokens. Discdirac / Deltaloop /
       Fixc_g / SpatialVec_cross all have this — and all four are REAL, so MakeNTKernel emits a
       single body and neither projection is ever built.
     * COMPLEX: `complexQ = !FreeQ[Diagrams, Complex]`, which makes MakeNTKernel emit all three
       bodies (Pure / RePart / complex) under `#if <FLOW>_VERDICT`. Every complex flow in tests/gen
       (ZA, ZA3, ZAqbq*, EtaPiL, Zq, ZAAqbq1Small, ...) is LINEAR in its trace tokens.

   With both, ntRePartIntegrand's linear decomposition `Coefficient[lin, tsym[t]]` returns a
   coefficient still containing the OTHER tokens' Unique["tr$"] placeholders — which are CForm'd into
   the C++ as bare `tr$2994` identifiers — and double-counts each multilinear term once per token it
   contains. ntPureIntegrand meanwhile drops the −Im(A)·Im(B) legs and COMPILES, i.e. is silently
   wrong. That is the bug this fixture pins. (Found via AqbqDirect tensors 5 and 6, whose projector
   puts its Lorentz index on a vec rather than a gamma; there the degree is 3.)

   THE INGREDIENTS, and why each is needed:
     anchor  a closed Dirac trace, tr(q1_slash ql_slash) = 4 (q1.ql)
     dirac2  a SECOND, index-disjoint closed Dirac trace — this is what disconnects the diagram
     lorT    a disjoint pure-Lorentz scalar, q3 . P_T(ql) . q3 — a third independent factor
     colT    tr(T^a T^b T^c) f^{abc}. The fundamental generators make the TRACE VALUES themselves
             complex (a real coefficient cannot: it would leave ntIm(tr) ≡ 0 and collapse RePart onto
             Pure). tr(TTT) = (d^{abc} + i f^{abc})/4 and d^{abc}f^{abc} = 0, so this contraction is
             PURELY IMAGINARY — which is why the overall factor of I below makes the integrand real
             and lands the probe on verdict 1 (RePart) rather than 0 (genuinely complex). This is the
             same shape the physics flows have; see the "all complex traces are pure-imaginary" note.
     I       sets complexQ, so all three bodies are emitted.

   The three diagrams give token degrees 1, 2 and 3 in one integrand, so the emitted RePart body
   contains BOTH the untouched linear path and the expansion — which is exactly the partition
   ntProjectIntegrand has to get right. (Measured: max degree 4 after the colour split, 2 of 3
   summands multilinear.)

   HOW IT IS GRADED. Nothing here is compared against a hand-derived closed form: the point is not
   the value but the AGREEMENT of the three emitted bodies. The untouched complex body is always
   correct, so it is the oracle, and compare_cplxdisc_num.cpp compiles this one kernel three times
   (once per verdict) and requires
       kernel<RePart>  ==  real(kernel<complex>)        exactly, to roundoff
   which is an IDENTITY, not a property of the flow — it must hold for any integrand, at any token
   degree. Reverting ntProjectIntegrand breaks it loudly (and also emits `tr$` into the header).
   The generated probe now asserts the same identity independently, over its own stubbed points.

   trailing '+' keeps the statement open across the line break (a leading '+' would end it) *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
                           Global`q1, Global`q2, Global`q3, Global`ql];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p};

anchor = ntGamma[m1, s1, s2] ntVec[Global`q1, m1] ntGamma[m2, s2, s1] ntVec[Global`ql, m2];
dirac2 = ntGamma[b1, e1, e2] ntVec[Global`q2, b1] ntGamma[b2, e2, e1] ntVec[Global`q3, b2];
lorT   = ntVec[Global`q3, i1] ntTransProj[Global`ql, i1, i2] ntVec[Global`q3, i2];
colT   = ntSUNT[3, a1, c1, c2] ntSUNT[3, a2, c2, c3] ntSUNT[3, a3, c3, c1] ntSUNf[3, a1, a2, a3];

net = I*anchor*colT +
      I*anchor*dirac2*colT +
      I*anchor*dirac2*lorT*colT;

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];

(* A real coefficient anywhere here would silently defeat the fixture: complexQ would go False, the
   three bodies would collapse to one, and the gate would pass while testing nothing. Assert it. *)
If[FreeQ[ntk[[1]]["Diagrams"], Complex],
  Print["FATAL: the fixture is not complex — MakeNTKernel would emit a single body and this gate ",
        "would test nothing. Check that the I and the fundamental-generator colour loop survived."];
  Exit[1]];

Print["[gen] Numeric backend ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_cplxdisc_num.cpp"}],
  FileNameJoin[{$gen, "Cplxdisc_num_kernel.hh"}],
  FileNameJoin[{$gen, "Cplxdisc_num_kernels.hh"}],
  "Name" -> "Cplxdisc_num_kernel", "Namespace" -> "cplxdisc_num",
  "Dressings" -> {}];

(* The whole point is that all THREE bodies were emitted and none of them carries a leaked scoped
   symbol. ntExportCpp's $-pattern already refuses the latter, but assert here too so a failure names
   this fixture rather than a generic emission guard. *)
Module[{h = Import[FileNameJoin[{$gen, "Cplxdisc_num_kernel.hh"}], "Text"]},
  If[StringCount[h, "VERDICT"] < 2,
    Print["FATAL: fewer than two VERDICT branches in the emitted kernel — the flow was not treated ",
          "as complex, so neither real projection was built."]; Exit[1]];
  If[StringContainsQ[h, "tr$"],
    Print["FATAL: leaked Unique[\"tr$\"] placeholder(s) in the emitted kernel."]; Exit[1]]];

(* The SAME ntk in single precision ("ComputeType" -> "float"), graded against the double complex
   body by compare_cplxdisc_float.cpp. Complex traces, the ntRe/ntIm split and a multilinear body make
   this the fixture with the most emitted types and literals for a stray double to hide in. *)
Print["[gen] Numeric backend, single precision ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_cplxdiscf_num.cpp"}],
  FileNameJoin[{$gen, "Cplxdiscf_num_kernel.hh"}],
  FileNameJoin[{$gen, "Cplxdiscf_num_kernels.hh"}],
  "Name" -> "Cplxdiscf_num_kernel", "Namespace" -> "cplxdiscf_num",
  "Dressings" -> {}, "ComputeType" -> "float"];

Module[{h = Import[FileNameJoin[{$gen, "Cplxdiscf_num_kernel.hh"}], "Text"],
        t = Import[FileNameJoin[{$gen, "Cplxdiscf_num_kernels.hh"}], "Text"]},
  If[StringContainsQ[h <> t, "double"],
    Print["FATAL: the float kernel still names `double` -- that arithmetic runs in double."]; Exit[1]]];

Print["Complex x multi-trace codegen kernels generated."];
