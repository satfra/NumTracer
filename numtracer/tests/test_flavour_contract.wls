#!/usr/bin/env wolframscript
(* Unit gate for the FUNDAMENTAL-FLAVOUR sector of the FunKit -> NumTracer translation:
   contractFlavour and promoteFlavResidue in mathematica/FunKitAdapter.m, and the leak guards
   that backstop them.

   WHY IT IS A UNIT TEST AND NOT ONLY A FLOW FIXTURE. flow_kaon_tadpole grades the emitted C++,
   but it lives behind NT_RUN_CODEGEN (minutes: Wolfram trace + generator compile + run) and so is
   effectively never run while editing FunKitAdapter.m. Everything below is pure symbolic algebra
   and takes about a second, so it runs in a PLAIN ctest.

   THE BUG THIS PINS. A flavour label rides along a line as a Kronecker delta. NumTracer has a cheap
   path (stitch chains, close a loop to Nf) and the SU(N) engine. It picked the engine only when it
   also saw flavour GENERATORS (TFlav), i.e. a pion/sigma Yukawa. A kaon-like field carries a
   fundamental flavour index but couples through no generator, so it got the cheap path — and in the
   tadpole, where the quartic gives delta_ac delta_bd + delta_ad delta_bc contracted against the
   loop's own delta, the deltas form a WEB, not a chain. Two rewrite rules cannot collapse a web, so
   the residue was CForm'd into the kernel as NumTracer_Private_flavDelta(F1, F2) — an undeclared
   identifier that both GCC and Clang parse.

   WHAT IS PINNED, and what breaks if it regresses:

     1. contractFlavour closes a chain in ALL FOUR index orientations. A Kronecker delta is
        symmetric, but the rule set was head-to-tail only, so a chain FunKit happened to emit in the
        other order silently failed to close. This is the "5 -> 4 unresolved deltas" the bug report
        measured for its candidate fix 2.

     2. delta^n folds to Nf for every n >= 2, not just n == 2. A surviving cube would be promoted
        and then read by labelCensus as a 3x-repeated label -> NumTrace::badlabel.

     3. THE BUG ITSELF: a delta WEB leaves no flavDelta behind, because promoteFlavResidue hands
        the residue to the SU(N) engine. Checked on the exact tadpole structure from
        repro_numtracer.wl.

     4. promoteFlavResidue is a strict NO-OP when the blind rules already close everything. This is
        what keeps the ~15 committed quark flows byte-identical, so it is the check that guards the
        whole design rather than the feature.

     5. The two flavour counts can never disagree. The blind path folds a loop to Global`Nf, the
        engine to the "FlavourGroup" rank; one diagram can use both, so a mismatch would silently mix
        conventions in a single coefficient (a wrong number, not a crash) — most likely in exactly an
        Nf = 2+1 setup whose light group is SU(2). FromFunKit must refuse.

     6. The emission backstop sees a private-context symbol. This is the class-level guard: any
        Mathematica private symbol reaching generated C++ is a leak, whichever package it escaped.

   Run standalone:  math -noprompt -script numtracer/tests/test_flavour_contract.wls  *)

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

Begin["NumTracer`Private`"];

(* Total number of check[] calls below. An Abort escaping mid-file would otherwise leave a passing
   tally simply because the remaining checks never executed; asserting the count closes that. *)
$nExpected = 26;
$fail = 0; nRan = 0;
check[label_, cond_] := (
  nRan++;
  Print[If[TrueQ[cond], "  PASS  ", "  FAIL  "], label];
  If[! TrueQ[cond], $fail++]);

SetAttributes[abortsQ, HoldRest];
abortsQ[expr_Hold, msg_] := CheckAbort[Quiet[ReleaseHold[expr]; False, msg], True];

(* The flavour heads are FormTracer's and live in Global` (verified with Names["*`deltaFundFlav"]
   before any reference). FromFunKit dispatches on the head NAME, so plain Global` symbols are a
   faithful stand-in for a TensorBases-produced expression — the same idiom the hand-built DSL
   fixtures in tests/gen use. *)
SetAttributes[{Global`deltaFundFlav, Global`TFlav}, {}];

Print["=== 1. contractFlavour: a Kronecker delta is SYMMETRIC ==="];
(* Every orientation of the shared index must close. Rules 2 and 3 are the ones that were missing. *)
check["head-to-tail   d[x,y] d[y,z] -> d[x,z]",
  contractFlavour[flavDelta[x, y] flavDelta[y, z]] === flavDelta[x, z]];
check["shared FIRST   d[y,x] d[y,z] -> d[x,z]",
  contractFlavour[flavDelta[y, x] flavDelta[y, z]] === flavDelta[x, z]];
check["shared LAST    d[x,y] d[z,y] -> d[x,z]",
  contractFlavour[flavDelta[x, y] flavDelta[z, y]] === flavDelta[x, z]];
check["closed loop    d[x,y] d[y,x] -> Nf",
  contractFlavour[flavDelta[x, y] flavDelta[y, x]] === Global`Nf];
check["self-trace     d[x,x]        -> Nf",
  contractFlavour[flavDelta[x, x]] === Global`Nf];
(* A three-link chain closing on itself, written in a deliberately mixed orientation. *)
check["mixed 3-chain  d[a,b] d[c,b] d[c,a] -> Nf",
  contractFlavour[flavDelta[a, b] flavDelta[c, b] flavDelta[c, a]] === Global`Nf];

Print[""];
Print["=== 2. contractFlavour: delta^n = delta, so it folds to Nf for every n >= 2 ==="];
check["d[x,y]^2 -> Nf", contractFlavour[flavDelta[x, y]^2] === Global`Nf];
check["d[x,y]^3 -> Nf", contractFlavour[flavDelta[x, y]^3] === Global`Nf];
check["d[x,y]^5 -> Nf", contractFlavour[flavDelta[x, y]^5] === Global`Nf];
(* A LONE open delta is not closeable and must survive for the engine — folding it would be wrong. *)
check["a lone open d[x,y] is NOT folded",
  contractFlavour[flavDelta[x, y]] === flavDelta[x, y]];

Print[""];
Print["=== 3. promoteFlavResidue: no-op when the blind rules already close everything ==="];
(* THE design check. If this ever fails, every committed quark flow's kernel moves. The identity
   must be on the LIST, not merely on its product: analyseDiagram numbers axis ids from the order of
   the tensor factors, so re-splicing would renumber them and change emitted code. *)
Block[{$ntFlavRank = 2},
  With[{fs = {ntGamma[v1, d1, d2], flavDelta[f1, f2], flavDelta[f2, f1], ntMetric[v1, v2]}},
    check["closed loop: factor list returned IDENTICALLY (===, not just same value)",
      promoteFlavResidue[fs] === fs]];
  With[{fs = {ntGamma[v1, d1, d2], ntMetric[v1, v2]}},
    check["no flavour delta at all: untouched", promoteFlavResidue[fs] === fs]];
  (* An OPEN delta cannot be closed by any rewrite, so it must be promoted into the engine. *)
  With[{fs = {ntGamma[v1, d1, d2], flavDelta[f1, f2]}},
    check["open delta is promoted to ntSUNDeltaFund[2,...]",
      FreeQ[promoteFlavResidue[fs], flavDelta] &&
      ! FreeQ[promoteFlavResidue[fs], ntSUNDeltaFund[2, f1, f2]]]]];

Print[""];
Print["=== 4. the bug: a delta WEB (the boson tadpole) leaves nothing uncontracted ==="];
(* The exact flavour structure of repro_numtracer.wl: the projector delta closing the external legs,
   the quartic's TWO pairings (an eager, purely SCALAR Plus that nothing distributes), and the
   tadpole loop's own delta. No chain rule can reach across that Plus, in any orientation. *)
web = flavDelta[Global`F1, Global`F2] *
      (flavDelta[Global`F1, Global`F45] flavDelta[Global`F49, Global`F2] +
       flavDelta[Global`F1, Global`F49] flavDelta[Global`F45, Global`F2]) *
      flavDelta[Global`F49, Global`F45];
check["the web genuinely defeats contractFlavour (else this test proves nothing)",
  ! FreeQ[contractFlavour[web], flavDelta]];
Block[{$ntFlavRank = 2},
  check["promoteFlavResidue leaves NO flavDelta behind",
    FreeQ[promoteFlavResidue[{web}], flavDelta]];
  check["...and routes it to the SU(2) engine",
    ! FreeQ[promoteFlavResidue[{web}], ntSUNDeltaFund[2, __]]]];
(* Unreachable in practice (flavDelta comes from exactly one $ffMap entry), but the promotion must
   never emit ntSUNDeltaFund[Automatic, ...]. *)
Block[{$ntFlavRank = Automatic},
  check["no rank available -> refuses instead of emitting a symbolic rank",
    abortsQ[Hold[promoteFlavResidue[{web}]], NumTrace::flavrank]]];

Print[""];
Print["=== 5. FromFunKit routes the whole fundamental-flavour sector, generators or not ==="];
Block[{Global`Nc = 3, Global`Nf = 2},
  (* WITHOUT any generator — the kaon case. This is what used to leak.
     Note the residue is EXPECTED to survive FromFunKit: contractFlavour runs twice, and the second
     pass (with the promotion) is in analyseDiagram, because a flavour chain buried in an eager
     dressed numerator only becomes a flat product after rewriteDressedNums. Pinning both halves
     documents that split — if the first check ever starts passing on its own, the promotion has
     been moved earlier and every dressed quark flow's kernel will have moved with it. *)
  With[{fk = FromFunKit[Global`FTerm[web /. flavDelta -> Global`deltaFundFlav]]},
    check["stage 1 (FromFunKit): the web's residue survives, as designed",
      ! FreeQ[fk, flavDelta]];
    Block[{$ntFlavRank = 2},
      check["stage 2 (promoteFlavResidue): end-to-end, no flavDelta is left",
        FreeQ[promoteFlavResidue[If[Head[fk] === Times, List @@ fk, {fk}]], flavDelta]]]];
  (* The control: the same algebra spelled through the generator route, which always worked.
     TFlav[0,f1,f2] = deltaFundFlav[f1,f2]/Sqrt[2 Nf], so the two spellings are the SAME object. *)
  check["the Sqrt[2 Nf] TFlav[0,..] spelling is clean too (control)",
    FreeQ[FromFunKit[Global`FTerm[web /. flavDelta[i_, j_] :>
            Sqrt[2 Global`Nf] Global`TFlav[0, i, j]]], flavDelta]];
  (* A closed flavour loop must still fold to the cheap scalar Nf = 2 and never reach the engine —
     this is the byte-identity of the committed quark flows, at FromFunKit level. *)
  check["a closed loop still folds cheaply to Nf (no SU(N) net built)",
    FreeQ[FromFunKit[Global`FTerm[
      Global`deltaFundFlav[f1, f2] Global`deltaFundFlav[f2, f1]]], ntSUNDeltaFund]]];

Print[""];
Print["=== 6. one flavour count, never two ==="];
(* The blind path folds a loop to Global`Nf; the engine folds it to the "FlavourGroup" rank. One
   diagram can use both routes, so if they disagree the coefficient mixes conventions silently. *)
Block[{Global`Nc = 3, Global`Nf = 3},
  check["FlavourGroup -> 2 while Global`Nf = 3 is refused",
    abortsQ[Hold[FromFunKit[Global`FTerm[Global`deltaFundFlav[f1, f2] Global`deltaFundFlav[f2, f1]],
                            "FlavourGroup" -> 2]], FromFunKit::flavcount]]];
Block[{Global`Nc = 3, Global`Nf = Global`NfSymbolic},
  check["unset Global`Nf (engine would silently default to 2) is refused",
    abortsQ[Hold[FromFunKit[Global`FTerm[Global`deltaFundFlav[f1, f2] Global`deltaFundFlav[f2, f1]]]],
            FromFunKit::flavcount]]];

Print[""];
Print["=== 7. the emission backstop sees a private-context symbol ==="];
(* CForm renders context marks as underscores, so a leaked private head prints as
   NumTracer_Private_flavDelta(F1, F2). It carries no List(...) tail and no $nnn, so the two generic
   patterns miss it — this is why it reached a kernel at all. *)
check["NumTracer_Private_flavDelta(F1, F2) is caught",
  AnyTrue[$ntCppLeakPatterns,
    StringPosition["  const auto _interp1 = NumTracer_Private_flavDelta(F1, F2);", #, 1] =!= {} &]];
check["the class generalises (FunKit_Private_foo)",
  AnyTrue[$ntCppLeakPatterns,
    StringPosition["  const double x = FunKit_Private_foo(a);", #, 1] =!= {} &]];
(* Ordinary emitted C++ must not trip it: no false positive on identifiers merely containing the
   words, nor on the runtime's own names. *)
check["no false positive on ordinary emitted C++",
  ! AnyTrue[$ntCppLeakPatterns,
    StringPosition["  const auto _interp1 = ZA(l1); int Private_count = 0; auto myPrivate_x = 1;",
                   #, 1] =!= {} &]];

Print[""];
If[nRan =!= $nExpected,
  Print["  FAIL  check count is ", nRan, ", expected ", $nExpected,
        " — a check was added/removed (update $nExpected) or an Abort escaped mid-file"];
  $fail++,
  Print["  PASS  all ", $nExpected, " checks ran (no early Abort)"]];
Print[""];
If[TrueQ[$fail > 0],
  Print["FAILED: ", $fail, " of ", nRan, " check(s)"]; End[]; Exit[1],
  Print["all ", nRan, " checks passed"]; End[]; Exit[0]];
