#!/usr/bin/env wolframscript
(* k-only lookup hoisting: ntPatchHoistWrappers (DiFfRG_compat.m) rewrites the wrapper TUs DiFfRG's
   MakeKernel scaffold writes, so they pass the M host-evaluated lookups the kernel now takes as extra
   parameters. The scaffold writes CT_get.cc as

       integrator.get(dest, <coordinates...>, k, <rest>);

   with ANY number of coordinates of ANY name (DiFfRG's "CoordinateArguments"), and two spaces before
   k. A wrapper left unpatched still compiles here; the mismatch only shows up when the MODEL is
   compiled, inside a DiFfRG integrator header ("expected 74, have 72"). So this pins (A) every
   coordinate shape the scaffold produces, and (B) that an unpatchable wrapper aborts at generation.

   Runs in a PLAIN ctest (no generation, no C++). *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

patch = NumTracer`Private`ntPatchHoistWrappers;
fail = False;
check[label_, got_, want_] := If[got === want,
  Print["  ok   ", label],
  (Print["FAIL: ", label, "\n  got:  ", got, "\n  want: ", want]; fail = True)];

(* One scaffold directory holding a CT_get.cc with the given call and the CT_map_1.cc every scaffold
   writes. Returns the patched CT_get.cc text, or "ABORTED". *)
run[getCall_String] := Module[{dir = CreateDirectory[], out},
  CreateDirectory[FileNameJoin[{dir, "src"}]];
  Export[FileNameJoin[{dir, "src", "CT_get.cc"}], "void f() {\n  " <> getCall <> "\n}\n", "Text"];
  Export[FileNameJoin[{dir, "src", "CT_map_1.cc"}],
    "auto g() {\n  return integrator.map(dest, coordinates, k, T, Zq);\n}\n", "Text"];
  out = CheckAbort[
    Quiet[patch[dir, "flow", 2], {NumTracer`Private`ntPatchHoistWrappers::unpatched}];
    Import[FileNameJoin[{dir, "src", "CT_get.cc"}], "Text"],
    "ABORTED"];
  DeleteDirectory[dir, DeleteContents -> True];
  out];

want[coords_String] := "void f() {\n  const auto _nth = flow_kernel<Regulator>::ntHoisted(k, T, Zq);\n" <>
  "  integrator.get(dest, " <> coords <> "k, T, Zq, _nth[0], _nth[1]);\n}";

(* --- A. every coordinate shape the scaffold writes -------------------------------------------- *)
check["no coordinate",          run["integrator.get(dest,  k, T, Zq);"],        want[""]];
check["one coordinate, p",      run["integrator.get(dest, p,  k, T, Zq);"],     want["p, "]];
check["one coordinate, p0",     run["integrator.get(dest, p0,  k, T, Zq);"],    want["p0, "]];
check["two coordinates, p0 p",  run["integrator.get(dest, p0, p,  k, T, Zq);"], want["p0, p, "]];
check["named integrator",       run["integrator_AD.get(dest,  k, T, Zq);"],
  StringReplace[want[""], "integrator.get" -> "integrator_AD.get"]];

(* --- B. a wrapper whose call the patch cannot recognise must abort, not be skipped ------------- *)
check["unrecognised get call aborts", run["integrator.get(dest, p, q, T, Zq);"], "ABORTED"];

If[fail, Print["hoist-wrapper test FAILED"]; Exit[1]];
Print["hoist-wrapper test passed"];
