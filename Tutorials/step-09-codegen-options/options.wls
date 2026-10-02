#!/usr/bin/env wolframscript
(* step-09 — Codegen options: what each one changes in the emitted C++.
 *
 * `MakeNTKernel` has close to thirty options. Almost all of them are supposed to be SEMANTICS-PRESERVING: they
 * change how the kernel is spelled — how diagrams are grouped, how subexpressions are shared, what
 * namespace and runtime it targets — without changing the number it returns.
 *
 * So the honest way to document them is to emit the SAME network several ways and check that every
 * variant agrees to the last bit. That is what this script sets up and the C++ checker verifies.
 * The differences you are meant to look at are in the generated source, not the output.
 *
 * Needs a Wolfram kernel and FunKit. No TensorBases, no FORM, no DiFfRG.
 *)

If[Quiet[Needs["NumTracer`"]] === $Failed,
  Module[{dir = Environment["NUMTRACER_DIR"]},
    If[dir === $Failed, Print["step-09: set NUMTRACER_DIR"]; Exit[1]];
    Get[FileNameJoin[{dir, "mathematica", "NumTracer.m"}]]]];
Global`Nc = 3;

$out = Environment["NT_TUTORIAL_OUT"];
If[$out === $Failed, $out = Directory[]];

{p, l1, cos1} = {Global`p, Global`l1, Global`cos1};
{q1, ql}      = {Global`q1, Global`ql};
myZ           = Global`myZ;

frame = propFrame[p, l1, cos1, q1, ql];
args  = {l1, cos1, p};
dress = {myZ};

(* ---------------------------------------------------------------------------------------------
 * @snip begin: net
 * A network with THREE diagrams — a top-level Plus is a sum of diagrams — two of which share the
 * same dressing coefficient myZ(l1). That shared coefficient is what the emitter factors
 * to collect.
 *
 * Closed form:   myZ(l1) * p^2   +   (q1.ql) * p^2   =   myZ(l1) p^2 + p^3 l1 cos1
 * (the first two diagrams are P^T + P^L = delta, so they add to q1.q1 = p^2).
 *
 * NOTE THE PARENTHESES — this is a genuine .wls trap. A line that is already a complete expression
 * ENDS there, so
 *     net = termA
 *         + termB;
 * assigns only termA and silently discards the rest. The kernel that comes out is valid, compiles,
 * and computes the wrong network. Wrapping the sum in ( ) is the fix; asserting the expected
 * DIAGRAM COUNT below is what catches it if you forget.
 *)
net = ( myZ[l1] ntVec[q1, mu]  ntTransProj[ql, mu, nu]   ntVec[q1, nu]
      + myZ[l1] ntVec[q1, rho] ntLongProj[ql, rho, sig]  ntVec[q1, sig]
      + ntSP[q1, ql] ntMetric[aa, bb] ntVec[q1, aa] ntVec[q1, bb] );
(* @snip end: net *)

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> dress];
Module[{nd = Length[ntk[[1]]["Diagrams"]]},
  Print["step-09: diagrams = ", nd];
  If[nd =!= 3, Print["step-09 FAIL: expected 3 diagrams, got ", nd]; Exit[1]]];

emit[nm_, ns_, opts___] := (
  MakeNTKernel[ntk,
    FileNameJoin[{$out, "gen_" <> ns <> ".cpp"}],
    FileNameJoin[{$out, nm <> ".hh"}],
    FileNameJoin[{$out, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> dress, opts];
  Print["step-09: emitted ", nm]);

(* ---------------------------------------------------------------------------------------------
 * @snip begin: variants
 * (1) The DEFAULT. Diagrams are grouped by their dressing coefficient and colour is folded
 *     numerically (both unconditional), so FunKit's emitter factors the shared myZ(l1) across the
 *     two projector diagrams instead of repeating it.
 *)
emit["Opt_default_kernel", "opt_default"];

(* (2) "Constant" -> expr. A loop-INDEPENDENT term flat-added to the integral. It does NOT enter
 *     the integrand: it is emitted as a separate `constant(...)` function that the consumer's
 *     integrator adds after quadrature. This is how an RG counterterm is supplied — step-13.
 *     Note the value returned by kernel() is unchanged; only constant() differs.
 *)
emit["Opt_const_kernel", "opt_const", "Constant" -> myZ[p]];

(* (3) Retargeting. A generated kernel is self-contained by DEFAULT: it includes only
 *     numtracer/codegen/runtime.hpp and lives in namespace `numtracer_kernels`. A consumer with its
 *     own support API redirects all of it. Here we only rename the namespaces — swapping
 *     "RuntimeInclude" too would need that header to exist, which is exactly what step-15 does when
 *     emitting against DiFfRG.
 *)
emit["Opt_ns_kernel", "opt_ns",
     "KernelNamespace" -> "step09_kernels", "SupportNamespace" -> "numtracer"];
(* @snip end: variants *)

Print["step-09: wrote 3 variants of the same kernel to ", $out];
