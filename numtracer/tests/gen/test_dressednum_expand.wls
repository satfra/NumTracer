#!/usr/bin/env wolframscript
(* dressedNumDecompose must split a propagator numerator over its TENSOR structure only.
   (a) finite-T psdash shape: Zq-dressed γ_μ (q_μ − c e0_μ) — the γ sum must become separate slash
       options, or the whole numerator stops being collectible (the reason for the Expand);
   (b) a composite scalar dressing −(Zq[q] + Zq[k] RF/|q|) must stay ONE option coefficient. Expanding
       it splits every atom into its own dressing slot (za4_147 kernel 3.1 MB → 14.9 MB).
   Synthetic input, no FunKit needed. *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

dec = NumTracer`Private`dressedNumDecompose;
nDN = NumTracer`Private`ntDressedNum;
ng = NumTracer`ntGamma; nd = NumTracer`ntDeltaDirac; nv = NumTracer`ntVec;

fails = 0;
expect[name_, cond_] := (Print["  ", name, ": ", If[TrueQ[cond], "ok", fails++; "FAIL"]]);

(* options of the (single) ntDressedNum in a decomposition *)
optsOf[d_] := First @ Cases[d, nDN[o_, _, _] :> o, {0, Infinity}];
(* round trip: the decomposition times out back to the input sum *)
roundTrip[d_, sum_] := Simplify[Expand[d /. nDN[o_, _, _] :> Total[(#[[1]] * Which[#[[2]] === {"ident"}, nd[i, j],
     True, Total[(#[[1]] * ng[mu, i, j] nv[#[[2]], mu]) & /@ #[[2, 2]]]]) & /@ o]] - Expand[sum]] === 0;

dress = -(Zq[q] + Zq[k] RF[k^2, q^2] / Sqrt[q^2]);

Print["== zero-T propagator numerator: composite dressing stays one option =="];
p0 = dress * ng[mu, i, j] nv[q, mu] + Mq[q] nd[i, j];
d0 = dec[p0];
expect["collectible", d0 =!= $Failed];
expect["2 options (ident, slash)", Length[optsOf[d0]] === 2];
expect["slash coefficient is the whole composite dressing",
  MemberQ[optsOf[d0][[All, 1]], c_ /; Simplify[c - dress] === 0]];
expect["round trip", roundTrip[d0, p0]];

Print["== finite-T psdash shape: the γ sum still splits =="];
pT = dress * ng[mu, i, j] (nv[q, mu] - c0 nv[e0, mu]) + Mq[q] nd[i, j];
dT = dec[pT];
expect["collectible", dT =!= $Failed];
expect["3 options (ident, slash q, slash e0)", Length[optsOf[dT]] === 3];
expect["each slash coefficient is (const) × the whole composite dressing",
  AllTrue[Select[optsOf[dT], #[[2]] =!= {"ident"} &][[All, 1]], FreeQ[Simplify[#/dress], Zq | RF] &]];
expect["round trip", roundTrip[dT, pT]];

Print[If[fails == 0, "\nALL OK", "\n" <> ToString[fails] <> " FAILED"]];
Exit[If[fails == 0, 0, 1]];
