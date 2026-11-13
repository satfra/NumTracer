#!/usr/bin/env wolframscript
(* step-23 — Trace tables: fill once per run, evaluate at every scale.
 *
 * A trace depends on the loop momentum l, the external momentum p and the angles. It is homogeneous
 * in (l, p), so once its denominators are factored out it is a short polynomial
 *
 *     N = sum_n  c_n(angles) l^n p^(D-n)
 *
 * whose coefficients c_n depend on the angles only. "TraceTables" makes MakeNTKernel emit those
 * coefficients as a table fill, and a second kernel entry point that rebuilds each trace from them.
 * A consumer fills the table once per run at its angle nodes; afterwards every evaluation, at any
 * l and p, reads a few numbers instead of recomputing the trace.
 *
 * This script emits one network twice: with every trace tabulated ("TraceTables" -> True), and with
 * NumTracer deciding per trace (Automatic). check_trace_tables.cpp fills the table by hand and
 * checks it against the inline kernel. Needs a Wolfram kernel and FunKit; no TensorBases, no FORM,
 * no DiFfRG.
 *)

If[Quiet[Needs["NumTracer`"]] === $Failed,
  Module[{dir = Environment["NUMTRACER_DIR"]},
    If[dir === $Failed, Print["step-23: set NUMTRACER_DIR"]; Exit[1]];
    Get[FileNameJoin[{dir, "mathematica", "NumTracer.m"}]]]];
Global`Nc = 3;

$out = Environment["NT_TUTORIAL_OUT"];
If[$out === $Failed, $out = Directory[]];

(* ---------------------------------------------------------------------------------------------
 * @snip begin: frame
 * The propagator frame of step-06. The kernel arguments are {l1, cos1, p}: the radial loop variable
 * FIRST, then the angle, then the external coordinate. That order is what the table needs: the
 * first argument is the radial scale, and the external coordinates (by default the argument named
 * p, set with "CoordinateArgs") follow the integration variables.
 *)
{p, l1, cos1} = {Global`p, Global`l1, Global`cos1};
{qp, ql}      = {Global`qp, Global`ql};
myZ           = Global`myZ;
frame = propFrame[p, l1, cos1, qp, ql];
args  = {l1, cos1, p};
(* @snip end: frame *)

(* ---------------------------------------------------------------------------------------------
 * @snip begin: net
 * Two diagrams. The first has transverse projectors on the loop momentum l and on l - p, so it
 * carries the denominators 1/l^2 and 1/(l-p)^2. The second is a Dirac trace times a dressing.
 * Closed forms, with c = cos1:
 *
 *   p.P^T(l).P^T(l-p).p              =  l p^2 (1 - c^2) (l - c p) / (l^2 - 2 c l p + p^2)
 *   tr[l-slash p-slash (l-p)-slash p-slash]  =  4 l^2 p^2 (2 c^2 - 1)  -  4 l p^3 c
 *
 * Note the parentheses around the sum: see step-09 for the .wls trap they avoid.
 *)
net = ( ntVec[qp, mu] ntTransProj[ql, mu, nu] ntTransProj[ql - qp, nu, rho] ntVec[qp, rho]
      + myZ[l1] ntVec[ql, a] ntVec[qp, b] ntVec[ql - qp, c] ntVec[qp, e]
          ntGamma[a, d1, d2] ntGamma[b, d2, d3] ntGamma[c, d3, d4] ntGamma[e, d4, d1] );
(* @snip end: net *)

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {myZ}];
Module[{nd = Length[ntk[[1]]["Diagrams"]]},
  Print["step-23: diagrams = ", nd];
  If[nd =!= 2, Print["step-23 FAIL: expected 2 diagrams, got ", nd]; Exit[1]]];

(* ---------------------------------------------------------------------------------------------
 * @snip begin: emit
 * The only new option. With True every trace is tabulated, and generation stops with an error if
 * the flow cannot be (a complex frame, a finite-temperature flow, or no consistent choice of
 * scales). Automatic tabulates a trace only where evaluating it inline costs more than reading its
 * table entries, and otherwise keeps it inline without complaint. The generator log reports the
 * outcome:
 *     [tables] 2 of 2 traces tabulated, 4 entries per table node
 *)
emit[nm_, ns_, tables_] := (
  MakeNTKernel[ntk,
    FileNameJoin[{$out, "gen_" <> ns <> ".cpp"}],
    FileNameJoin[{$out, nm <> ".hh"}],
    FileNameJoin[{$out, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> {myZ},
    "TraceTables" -> tables];
  Print["step-23: emitted ", nm]);

emit["Tab_kernel", "tab_all", True];
emit["Tab_auto_kernel", "tab_auto", Automatic];
(* @snip end: emit *)
