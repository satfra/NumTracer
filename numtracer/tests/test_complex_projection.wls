#!/usr/bin/env wolframscript
(* Unit gate for the real/imaginary projections of a complex integrand — ntProjectIntegrand and
   friends in mathematica/CodegenRealProjection.m.

   WHY IT IS A UNIT TEST AND NOT A FLOW FIXTURE. Every other check of this machinery (cplxdisc_num,
   and the numeric probe that every complex flow runs) goes through a full generation: Wolfram
   tracing, a C++ generator, a compile. That takes minutes, so it lives behind NT_RUN_CODEGEN and is
   effectively never run while editing Codegen*.m. The properties below are pure symbolic algebra and
   take under a second, so this runs in a PLAIN ctest — it is the first fast test that can see a
   change to mathematica/.

   It complements rather than replaces cplxdisc_num: that fixture grades the emitted C++
   (kernel<RePart> == real(kernel<complex>)); this one grades the symbolic step that produces it,
   including the two failure modes that never reach C++ at all because they abort first.

   WHAT IS PINNED, and what breaks if it regresses:

     1. ntImagSplitSafeQ — which coefficient shapes the I -> ii substitution may be applied to. This is a
        SILENT-WRONG-ANSWER guard, not a style check: Coefficient power-series expands anything that
        is not a polynomial in ii and returns the leading term, so "tr0"/(a + I b) projects to
        ntRe["tr0"]/a. At finite density that deletes muq from every quark propagator and the kernel
        compiles clean. If this table loosens, that comes back.

     2. The one-pass ntRePartLinear reproduces the per-token-Coefficient formulation EXACTLY (===),
        not merely up to Simplify. That exactness is what keeps every committed kernel
        byte-identical, so the reference implementation is inlined below and compared structurally.

     3. Both projections — RePart AND Pure — refuse a nested Complex. ntPureLinear does not go
        through ntSplitRealImag and has its own copy of the blind spot, so a guard placed inside the
        split would cover only half of them.

     4. The multilinear (disconnected-diagram) route stays exact against ComplexExpand at degree 3.

   Run standalone:  math -noprompt -script numtracer/tests/test_complex_projection.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

(* Total number of check[] calls below. An Abort escaping mid-file would otherwise leave a passing
   tally simply because the remaining checks never executed; asserting the count closes that. *)
$nExpected = 34;
$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);

(* Abort[] is how every guard in Codegen*.m reports; CheckAbort turns that into a value so one
   refusal does not tear down the rest of the file. Quiet suppresses only the message we expect.
   HoldRest is load-bearing: without it `MakeNTKernel::cplxnest` evaluates to the message TEXT before
   Quiet ever sees it, Quiet rejects it (Quiet::mgre) and the abort is not caught — the check then
   reports FAIL against working code. *)
SetAttributes[abortsQ, HoldRest];
abortsQ[expr_Hold, msg_] :=
  CheckAbort[Quiet[ReleaseHold[expr]; False, msg], True];

Print["=== 1. ntImagSplitSafeQ: which shapes the ii-substitution may be applied to ==="];
(* SAFE: a Complex itself, and anything built from Complexes by Times / Plus / non-negative integer
   Power — exactly the shapes Coefficient[·, ii, n] handles as a genuine polynomial. *)
Scan[check["safe:   " <> ToString[InputForm[#]], ntImagSplitSafeQ[#]] &,
  {x*y, I*x*y, Complex[3, -2]*x, (1 + I a)*(2 + I b), I*(x + I y), (a + I b)^3,
   Zq/(l0 + Ep), 1/(a + b)^2}];
(* UNSAFE: a Complex under a negative power (every denominator — the finite-density l0 + I muq) or
   under a non-rational head, where Coefficient silently returns a series leading term. *)
Scan[check["unsafe: " <> ToString[InputForm[#]], ! ntImagSplitSafeQ[#]] &,
  {1/(a + I b), 1/(a + I b)^2, Zq/(l0 + I muq + Ep), Sqrt[a + I b], Exp[I a],
   1/(a + I b) + 1/(c - I d)}];
(* the message has to name the DENOMINATOR, not the whole summand it is buried in *)
check["ntImagSplitUnsafeParts isolates the offending factor",
  ntImagSplitUnsafeParts[Zq*dr1*"x"/(l0 + I muq + Ep)] === {1/(Ep + l0 + I muq)}];

Print[""];
Print["=== 2. one-pass ntRePartLinear === the per-token-Coefficient reference ==="];
(* The implementation this replaced, verbatim. Kept HERE rather than in Codegen*.m: its only remaining
   purpose is to be the oracle for the byte-identity claim. *)
refRePartLinear[integrand_] := Module[{toks = Union[Cases[integrand, _String, Infinity]], tsym, lin, ii},
  tsym = AssociationThread[toks -> Table[Unique["tr$"], {Length[toks]}]];
  lin = integrand /. tsym;
  Total[Function[t, Module[{cc = Coefficient[lin, tsym[t]] /. Complex[ar_, ai_] :> ar + ii*ai},
     Coefficient[cc, ii, 0]*Global`ntRe[t] - Coefficient[cc, ii, 1]*Global`ntIm[t]]] /@ toks]
    + ((lin /. Thread[Values[tsym] -> 0]) /. Complex[ar_, ai_] :> ar)];

SeedRandom[42];
dress = Table[Symbol["dr" <> ToString[i]], {i, 12}];
den = Table[Symbol["dn" <> ToString[i]], {i, 6}];
mkSummand[nT_] :=
  RandomChoice[{1, -1, 2, I, -I, Complex[3, -2]}]*(Times @@ RandomSample[dress, 4])*
    RandomChoice[den]^RandomChoice[{-1, -2, -3}]*RandomInteger[{1, 999}]*
    ("tr" <> ToString[RandomInteger[{0, nT - 1}]]);
(* token-free summands (real and imaginary) exercise the constant bucket *)
Scan[
  Function[nT,
    Module[{ig = Total[Table[mkSummand[nT], {600}]] + dr1*dn1 + I*dr2},
      check["structurally identical at " <> ToString[nT] <> " tokens",
        ntRePartLinear[ig] === refRePartLinear[ig]]]],
  {1, 5, 20, 77}];
(* hand-checked small case, so a co-regression of both implementations cannot pass *)
check["Re(I Zq tr0 + 3 tr1) == -Zq ntIm(tr0) + 3 ntRe(tr1)",
  ntRePartLinear[I*Zq*"tr0" + 3*"tr1"] === -Zq*Global`ntIm["tr0"] + 3*Global`ntRe["tr1"]];
check["Pure drops the imaginary coefficient",
  ntPureIntegrand[I*Zq*"tr0" + 3*"tr1"] === 3*Global`ntRe["tr1"]];

Print[""];
Print["=== 3. nested Complex is refused by BOTH projections ==="];
badDen = "tr0"*Zq/(l0 + I muq + Ep);
check["ntRePartIntegrand aborts (cplxnest)",
  abortsQ[Hold[ntRePartIntegrand[badDen]], MakeNTKernel::cplxnest]];
check["ntPureIntegrand aborts (cplxnest)",
  abortsQ[Hold[ntPureIntegrand[badDen]], MakeNTKernel::cplxnest]];
(* no false positive: a REAL denominator with a complex numerator is the ordinary flow shape *)
check["real denominator still projects",
  ntRePartIntegrand[I*Zq*"tr0"/(l0 + Ep) + 3*"tr1"] ===
    -((Zq*Global`ntIm["tr0"])/(Ep + l0)) + 3*Global`ntRe["tr1"]];
(* a token that is not a bare factor was silently DROPPED by the old Coefficient extraction *)
check["nested trace token aborts (toknest)",
  abortsQ[Hold[ntRePartLinear[c*"tr0"^2]], MakeNTKernel::toknest]];

Print[""];
Print["=== 4. multilinear route stays exact at degree 3 ==="];
mig = I*Zq*"trA"*"trB"*"trC" + 2*"trA" + Complex[3, -2]*"trB"*"trC" + 7;
sub = {Global`ntRe["trA"] -> ra, Global`ntIm["trA"] -> ia, Global`ntRe["trB"] -> rb,
       Global`ntIm["trB"] -> ib, Global`ntRe["trC"] -> rc, Global`ntIm["trC"] -> ic};
check["RePart == ComplexExpand[Re[...]]",
  Expand[(ntRePartIntegrand[mig] /. sub) -
    ComplexExpand[Re[mig /. {"trA" -> ra + I ia, "trB" -> rb + I ib, "trC" -> rc + I ic}]]] === 0];
(* the placeholders that leaked into emitted C++ as bare `tr$2994` identifiers *)
check["no scoped symbol survives either route",
  Cases[{ntRePartIntegrand[mig], ntPureIntegrand[mig], ntRePartLinear[I*Zq*"tr0" + 3*"tr1"]},
        s_Symbol /; StringContainsQ[SymbolName[s], "$"], {0, Infinity}] === {}];

Print[""];
Print["=== 5. \"ComplexRuntimeProjection\": finite-density denominators ==="];
check["option declared on MakeNTKernel",
  MemberQ[Keys[Association @@ Options[NumTracer`MakeNTKernel]], "ComplexRuntimeProjection"]];
check["option declared on mkGenerateKernel (so MakeNTKernel's FilterRules forwards it)",
  MemberQ[Keys[Association @@ Options[mkGenerateKernel]], "ComplexRuntimeProjection"]];
$ntComplexRuntimeProjection = True;
den = Zq/(l0 + I muq + Ep);
(* Re(z t) = Re(z)Re(t) - Im(z)Im(t), with z left INTACT for the generated code to split. Keeping it
   factored is the whole point — expanding it is what defeated CSE and made generation stall. *)
check["RePart defers the split to runtime, denominator unexpanded",
  ntRePartIntegrand[den*"tr0" + 2*"tr1"] ===
    -(Global`ntIm["tr0"]*Global`ntIm[den]) + 2*Global`ntRe["tr1"] +
      Global`ntRe["tr0"]*Global`ntRe[den]];
check["Pure drops the -Im*Im leg",
  ntPureIntegrand[den*"tr0" + 2*"tr1"] ===
    2*Global`ntRe["tr1"] + Global`ntRe["tr0"]*Global`ntRe[den]];
(* Numerical proof that the deferred form is the true real part — evaluate ntRe/ntIm as the actual
   complex operations and compare against Re of the untouched integrand. Guards against a sign or a
   swapped leg, which the structural checks above would accept if BOTH sides were written wrong. *)
Module[{vals = {Zq -> 1.7, l0 -> 0.3, muq -> 0.8, Ep -> 2.1}, tr0v = 0.4 - 1.3 I, exact, viaRt},
  exact = Re[(den*tr0v + 2*(0.9 + 0.2 I)) /. vals];
  viaRt = ntRePartIntegrand[den*"tr0" + 2*"tr1"] /.
      {Global`ntRe[s_String] :> Re[If[s === "tr0", tr0v, 0.9 + 0.2 I]],
       Global`ntIm[s_String] :> Im[If[s === "tr0", tr0v, 0.9 + 0.2 I]],
       Global`ntRe[z_] :> Re[z /. vals], Global`ntIm[z_] :> Im[z /. vals]};
  check["runtime form reproduces Re(integrand) numerically", Abs[viaRt - exact] < 1.*^-12]];
(* the OLD symbolic path on the same input, to show what the guard is protecting against *)
check["the symbolic split would have deleted muq (this is the bug)",
  Module[{ii}, Coefficient[den /. Complex[ar_, ai_] :> ar + ii*ai, ii, 0]] === Zq/(Ep + l0)];
$ntComplexRuntimeProjection = False;
check["flag restored: nested Complex aborts again",
  abortsQ[Hold[ntRePartIntegrand[den*"tr0"]], MakeNTKernel::cplxnest]];

(* The verdict MUST be taken inside the context $fail was assigned in. Taken after End[] it reads an
   undefined Global`$fail, `$fail > 0` evaluates to neither True nor False, If returns unevaluated,
   NEITHER branch runs and the script falls off the end with exit 0 — a green test that ran no
   verdict at all. *)
(* NOT routed through check[]: its argument is evaluated before check[] bumps the counter, so the
   expected number would have to be one less than the total it reports — an off-by-one waiting to be
   "fixed" wrongly by the next person. Asserted directly instead. *)
If[nRan =!= $nExpected,
  Print["  FAIL  check count is ", nRan, ", expected ", $nExpected,
        " — a check was added/removed (update $nExpected) or an Abort escaped mid-file"];
  $fail++,
  Print["  PASS  all ", $nExpected, " checks ran (no early Abort)"]];
Print[""];
If[TrueQ[$fail > 0],
  Print["FAILED: ", $fail, " of ", nRan, " check(s)"]; End[]; Exit[1],
  Print["all ", nRan, " checks passed"]; End[]; Exit[0]];
