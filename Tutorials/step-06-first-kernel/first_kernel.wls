#!/usr/bin/env wolframscript
(* step-06 — Your first generated kernel.
 *
 * This is the smallest possible end-to-end use of the Mathematica front-end. It generates a C++
 * kernel for exactly the network step-04 contracted by hand,
 *
 *     p . P(l) . p / p^2   =   1 - cos^2(theta),
 *
 * and writes it to disk as an ordinary, self-contained header. The network is hand-built, so
 * TensorBases, FORM and DiFfRG are not involved; MakeNTKernel itself writes the C++ through
 * FunKit's emitter, so FunKit must be installed.
 *
 * Run it either through the tutorial CMake project
 *     cmake --build build --target first_kernel_gen
 * or directly
 *     NT_TUTORIAL_OUT=/tmp/out wolframscript -script first_kernel.wls
 *)

(* ---------------------------------------------------------------------------------------------
 * @snip begin: load
 * Locating the front-end. When NumTracer is installed (`cmake --install`), `Needs["NumTracer`"]`
 * just works because the installer puts the paclet on Wolfram's application path. The NUMTRACER_DIR
 * fallback is what lets this script also run straight from a repository checkout — the tutorial
 * CMake project sets it.
 *)
If[Quiet[Needs["NumTracer`"]] === $Failed,
  Module[{dir = Environment["NUMTRACER_DIR"]},
    If[dir === $Failed,
      Print["step-06: NumTracer not found. Install it, or set NUMTRACER_DIR to the numtracer/ ",
            "directory of a checkout."];
      Exit[1]];
    Get[FileNameJoin[{dir, "mathematica", "NumTracer.m"}]]]];

(* NumTrace needs a colour rank even for a pure-Lorentz network: every SU(N) head carries its rank
 * as a leading argument, and the analyser checks them up front. Nothing here uses colour. *)
Global`Nc = 3;
(* @snip end: load *)

$out = Environment["NT_TUTORIAL_OUT"];
If[$out === $Failed, $out = Directory[]];

(* ---------------------------------------------------------------------------------------------
 * @snip begin: frame
 * The frame. `propFrame[p, l1, cos1, qp, ql]` declares the one-angle propagator kinematics:
 * the momentum tagged `qp` lies along axis 0 with magnitude `p`, and the loop momentum tagged
 * `ql` has magnitude `l1` at angle `cos1` to it. Those three symbols — p, l1, cos1 — are the
 * kernel's runtime arguments.
 *
 * Note the two kinds of name. `qp`/`ql` are MOMENTUM TAGS: they appear inside the network and are
 * resolved to components by the frame. `p`/`l1`/`cos1` are SCALARS: they appear in the generated
 * C++ signature. Keeping them distinct is the single most common source of confusion in a first
 * script (see step-07).
 *)
{p, l1, cos1} = {Global`p, Global`l1, Global`cos1};   (* runtime scalar kernel arguments *)
{qp, ql}      = {Global`qp, Global`ql};               (* momentum tags: external, loop *)

frame = propFrame[p, l1, cos1, qp, ql];
args  = {l1, cos1, p};
(* @snip end: frame *)

(* ---------------------------------------------------------------------------------------------
 * @snip begin: net
 * The network, in the DSL. Repeated index labels (mu, nu) are summed — Einstein convention, the
 * same rule as the C++ side. `ntSP[qp, qp]` is the scalar product qp.qp = p^2; it is a SCALAR
 * coefficient, not a tensor factor, so it multiplies the network rather than contracting with it.
 *
 * This is exactly the step-04 network, divided by p^2 so that the answer is the pure angular
 * factor 1 - cos^2(theta).
 *)
net = (1/ntSP[qp, qp]) ntVec[qp, mu] ntTransProj[ql, mu, nu] ntVec[qp, nu];
(* @snip end: net *)

(* ---------------------------------------------------------------------------------------------
 * @snip begin: emit
 * The two entry points.
 *
 *   NumTrace[net, "Frame" -> …, "Args" -> …]   analyses the network into an NTKernel: the list of
 *     diagrams (each a scalar coefficient times a contraction), the env-id layout, and the frame.
 *     Nothing is contracted yet — this is analysis only.
 *
 *   MakeNTKernel[ntk, genFile, kernelFile, tracesFile, …]   does the rest: it EMITS a small C++
 *     generator program to genFile, COMPILES it, RUNS it (which is where the numeric contraction
 *     of steps 01-05 actually happens), and writes the two committed headers.
 *
 * Three output files, always. The `traces` header holds the lowered `trN(const double* f)`
 * functions — the straight-line arithmetic step-05 printed. The `kernel` header holds the class
 * that fills the frame symbols and calls them.
 *)
ntk = NumTrace[net, "Frame" -> frame, "Args" -> args];

Print["step-06: diagrams = ", Length[ntk[[1]]["Diagrams"]]];

MakeNTKernel[ntk,
  FileNameJoin[{$out, "gen_first_kernel.cpp"}],   (* the generator program (a build artefact) *)
  FileNameJoin[{$out, "first_kernel.hh"}],        (* the kernel class *)
  FileNameJoin[{$out, "first_kernels.hh"}],       (* the lowered traces *)
  "Name"      -> "first_kernel",
  "Namespace" -> "step06",
  "Dressings" -> {}];
(* @snip end: emit *)

Print["step-06: wrote first_kernel.hh / first_kernels.hh to ", $out];
