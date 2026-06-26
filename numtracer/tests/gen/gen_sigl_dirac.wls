#!/usr/bin/env wolframscript
(* Regression test for the hSigL meson-sector δ-closed-loop collapse (BUG_meson_sector_dressed_collapse.md).

   The σL (scalar) component of the quark self-energy is projected by the spinor IDENTITY δ[d1,d2]
   (the σL leg has no Dirac tensor structure), so the qbar-q line is closed by a BARE spinor δ sitting
   INSIDE the Dirac chain. A meson-exchange contribution threads two pion γ5 vertices, giving the loop
       tr( δ[d1,d2] · γ5 · S · γ5 · S ) = tr(γ5·(γ·l1)·γ5·(γ·p1)) ∝ l1·p1  ≠ 0  (a real, kinematic trace).

   The old `orderDiracFacs` walked the spinor loop with a DIRECTED forward map (`byIn`); the symmetric
   external δ ran parallel to the chain, collided on the shared in/out keys and broke the walk, so the
   trace collapsed to a single γ5 → tr(γ5)=0 → the whole diagram vanished. The fixed walk is undirected
   and closes the cycle. This flow has NO other in-repo coverage of a δ-closed loop; it is the permanent
   guard. The Dense DTensor oracle contracts the network entry-for-entry (no ordering walk), so its
   value is unaffected by the bug — numeric (the gamma-trace theorem) must match it.

   Hand-built (no FunKit derivation): a closed σL-style δ loop expressed directly in the FunKit-token
   vocabulary and lowered through FromFunKit + NumTrace, exactly like a real traced flow.

   Run from numtracer/:   wolfram -script tests/gen/gen_sigl_dirac.wls *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

(* colour rank baked at import time; this flow carries no colour/flavour heads (pure Dirac), so the
   rank is irrelevant — set it for FromFunKit's sunMap and pass a trivial flavour group. *)
SetNc[3];

(* baked closed-form fermion dressing G(q^2) = 1/(q^2 + R_F(k^2,q^2)); RF supplied by the C++ Reg. *)
G[q_] := 1/(sp[q, q] + RF[k^2, sp[q, q]]);

(* the open quark line d1→a→b→c→d2: γ5 (pion), q-slash propagator (loop momentum l1), γ5 (pion),
   q-slash propagator (external momentum p1). Each propagator carries the usual (-I) γ·q numerator. *)
line = gamma5[d1, a] *
       ((-I) gamma[m1, a, b] vec[l1, m1]) * G[l1] *
       gamma5[b, c] *
       ((-I) gamma[m2, c, d2] vec[p1, m2]) * G[p1];

(* external σL scalar projector: the bare spinor δ[d1,d2] closing the line (coeff 1/4 = the qbar-q
   normalisation 1/(4·tr 1) stand-in; the value is arbitrary, the TOPOLOGY is the point). *)
myTrace = FEx[ FTerm[1/4, deltaDirac[d1, d2]] ** line ];

net   = FromFunKit[myTrace, "FlavourGroup" -> 2, "DressingCollection" -> True];
frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
ntk   = NumTrace[net, "Frame" -> frame, "Args" -> {Global`l1, Global`cos1, Global`p, Global`k}];

MakeNTKernel[ntk, FileNameJoin[{$gen, "Sigl_dirac_dense_kernel.hh"}],
  "Backend" -> "Dense", "Name" -> "Sigl_dirac_dense", "Namespace" -> "sigl_dirac_dense"];
MakeNTKernel[ntk, FileNameJoin[{$gen, "gen_sigl_dirac_num.cpp"}],
  FileNameJoin[{$gen, "Sigl_dirac_num_kernel.hh"}], FileNameJoin[{$gen, "Sigl_dirac_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Sigl_dirac_num_kernel", "Namespace" -> "sigl_dirac_num"];
Print["σL δ-closed meson-loop regression kernels generated (Dense oracle + numeric)."];
