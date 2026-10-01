(* ---- the two REAL projections of a complex integrand ---------------------------------------
   Both are emitted for every complex flow; the imaginary-part probe picks between them (and the
   untouched complex form) with a preprocessor `#if` — see ntProbeSource. Which one is valid is a
   numerical question about the generated traces, but BUILDING them is pure symbolic rewriting, so
   both happen here, before any C++ exists.

   "Pure": the projection Complex -> Re is exact, i.e. Σ Im(c)·tr ≈ 0. Wrap the trace tokens in ntRe
   FIRST so the kernel is provably double-typed even when a trace function is complex-typed, then drop
   the imaginary coefficients outright. This is the cheap body — the imaginary half never appears. *)
ntPureLinear[integrand_] := (integrand /. s_String :> Global`ntRe[s]) /. Complex[a_, b_] :> a;

(* ---- the token-degree problem, and the multilinear generalisation --------------------------
   BOTH real projections below were written against a stated precondition — "the integrand is LINEAR
   in the trace tokens" — that the assembly does not honour. A DISCONNECTED diagram (lorFacOf =!= None:
   a Dirac trace times one or more disconnected pure-Lorentz scalar components) contributes
   `Π traceRef[factor groups] · traceRef[anchor]`, i.e. a PRODUCT of trace tokens — see the integrand
   Sum near the end of mkGenerateKernel. Measured degree 3 on the AqbqDirect T5/T6 quark-gluon
   projections, whose Lorentz index sits on a vec rather than a gamma.

   What that broke, all three silently or loudly wrong:
     * ntRePartLinear's `Coefficient[lin, tsym[t]]` returns a coefficient that STILL CONTAINS the other
       tokens' `Unique["tr$"]` placeholders, which are then CForm'd into the C++ as bare `tr$2994`.
     * the same term is visited once per token it contains, so it is DOUBLE (triple, …) COUNTED.
     * a repeated token (t^2) is dropped outright: Coefficient[c t^2, t] == 0 and the token-free
       remainder kills it too.
     * ntPureLinear maps each token to ntRe independently, so `Re(c·A·B)` loses the −Im(A)·Im(B) leg.
       That one COMPILES, i.e. it is a wrong number with no diagnostic.

   The generalisation is per-SUMMAND and exact for any degree: split a summand's factors into the
   string-free part (left FACTORED — this is what keeps COEN's CSE alive, and why ComplexExpand is
   still not used anywhere here) and the token-bearing part, expand ONLY the latter (it is a product
   of two or three short sums), then take the real part of each resulting monomial with the same
   `ii`-substitution trick the linear path uses.

   COST: for degree n the exact real part is 2^n products against the broken code's 2n, i.e. IDENTICAL
   at n = 2 and +33% at n = 3, on the multilinear summands only. *)

ntSummandsOf[e_] := If[Head[e] === Plus, List @@ e, {e}];
ntFactorsOf[e_] := If[Head[e] === Times, List @@ e, {e}];

(* the trace tokens of one MONOMIAL, with multiplicity (t^2 counts twice) *)
ntTokensOfMonomial[m_] :=
  Flatten[
    Function[f,
        Which[
          StringQ[f], {f},
          MatchQ[f, Power[_String, _Integer?Positive]], ConstantArray[f[[1]], f[[2]]],
          True, {}]] /@ ntFactorsOf[m]];

(* Expand ONLY the string-bearing factors of a summand. Returns {stringFreePart, {monomials..}};
   the string-free part is never expanded. *)
ntSplitTokenPart[s_] := Module[{facs = ntFactorsOf[s], tokFacs, plain},
  tokFacs = Select[facs, ! FreeQ[#, _String] &];
  plain = Times @@ Select[facs, FreeQ[#, _String] &];
  If[tokFacs === {}, {plain, {}}, {plain, ntSummandsOf[Expand[Times @@ tokFacs]]}]];

(* highest number of trace tokens multiplied together anywhere in this summand *)
ntTokenDegree[s_] := Module[{monos = Last[ntSplitTokenPart[s]]},
  If[monos === {}, 0, Max[Length[ntTokensOfMonomial[#]] & /@ monos]]];

(* ii^2 -> -1 on an expression already expanded in the real stand-in `ii` for the imaginary unit *)
ntIiReduce[e_, ii_] := e //. Power[ii, n_Integer /; n >= 2] :> (-1)^Quotient[n, 2] * ii^Mod[n, 2];

(* {Re[e], Im[e]} for a string-free coefficient whose only imaginary content is explicit `Complex`
   numbers (which is all the assembly ever produces — no symbol here is imaginary).
   The two cheap shapes are taken WITHOUT expanding, and they cover essentially every real coefficient
   and every "one Complex scaling a factored real expression": that is what keeps the dressing
   coefficients factored for COEN's CSE. Anything else falls back to the same ii-substitution the
   linear path has always used, so the factoring behaviour there is no worse than before. *)
(* ---- when is the ii-substitution faithful? --------------------------------------------------
   Substituting I -> ii (a real symbol) and reading off Coefficient[·, ii, 0|1] is EXACT only where
   the expression is a POLYNOMIAL in ii. Mathematica does not complain otherwise — it power-series
   expands and hands back the leading term. Measured consequences on shapes the finite-density
   propagator actually produces:

     "tr0"/(a + I b)   ->  ntRe["tr0"]/a        the imaginary part of the DENOMINATOR is deleted
     "tr0"/(a + I b)^2 ->  ntRe["tr0"]/a^2      (same; at finite mu, muq vanishes from every quark
                                                 propagator and the kernel is a wrong number that
                                                 compiles — only the numeric probe can catch it)
     Sqrt[a + I b]     ->  Sqrt[a + b ii$7182]  the whole expression comes back as the ii^0 part,
                                                 so the Module-local stand-in LEAKS into the C++

   It is also where the cost is: with 1/(l0 + I muq + E)^n factors the projection went 1.5 s at 50
   summands, 15 s at 100, 39 s at 200, 62 s at 400 — the reported "generation stalls after complex
   projection". The per-token Coefficient loop this file used to run is a rounding error next to it
   (0.13 s at 77 tokens over 4000 summands).

   The test is purely STRUCTURAL — deliberately no Expand, so it costs nothing on a large factored
   coefficient, which is the whole point of keeping those factored. It admits exactly the shapes the
   substitution can handle: a Complex itself, products and sums of admissible parts, and a
   non-negative integer power of one. Everything else that carries a Complex is refused. *)
ntIiSafeQ[e_] :=
  FreeQ[e, Complex] ||
  Switch[Head[e],
    Complex, True,
    Times | Plus, AllTrue[List @@ e, ntIiSafeQ],
    Power, IntegerQ[e[[2]]] && NonNegative[e[[2]]] && ntIiSafeQ[e[[1]]],
    _, False];

(* The MAXIMAL offending nodes, for the message: descend through the Times/Plus skeleton (whose parts
   are individually admissible or not) and stop at the first node that is not. Reporting the whole
   summand instead would bury the one denominator that matters in a page of factored dressing. *)
ntIiUnsafeParts[e_] :=
  If[ntIiSafeQ[e], {},
    Switch[Head[e],
      Times | Plus, Flatten[ntIiUnsafeParts /@ (List @@ e)],
      _, {e}]];

ntSplitRealImag[e_] := Module[{fs, cs, rest, c, ii, t},
  If[FreeQ[e, Complex], Return[{e, 0}]];
  fs = ntFactorsOf[e];
  cs = Cases[fs, _Complex];
  rest = DeleteCases[fs, _Complex];
  If[cs =!= {} && FreeQ[rest, Complex],
    c = Times @@ cs;
    Return[{Re[c] * Times @@ rest, Im[c] * Times @@ rest}]];
(* Not a polynomial in the stand-in: the extraction below would silently return a series leading term
   (see ntIiSafeQ). Under "ComplexRuntimeProjection" hand the untouched, still-FACTORED expression to
   the generated code and let it take the parts at runtime. Otherwise ntProjectIntegrand has already
   refused, so this is unreachable and the extraction below keeps its historical behaviour. *)
  If[! ntIiSafeQ[e] && TrueQ[$ntComplexRuntimeProjection],
    Return[{Global`ntRe[e], Global`ntIm[e]}]];
  t = ntIiReduce[Expand[e /. Complex[ar_, ai_] :> ar + ii*ai], ii];
  {Coefficient[t, ii, 0], Coefficient[t, ii, 1]}];

(* {Re, Im} of Σ_m coef_m · Π_j t_j, each token replaced by ntRe[t] + i·ntIm[t]. Only this small
   token polynomial is expanded — degree 2-3, and its per-monomial coefficients are the factor-group
   dressing tokens (1 in every undressed flow). *)
ntTokenProductRealImag[monos_List, pureQ_] := Module[{ii, tot},
  tot =
    Total[
      Function[m,
          Module[{cr, ci},
            {cr, ci} = ntSplitRealImag[Times @@ Select[ntFactorsOf[m], FreeQ[#, _String] &]];
            If[TrueQ[pureQ], ci = 0];
            (cr + ii*ci) * Times @@ (Function[t, Global`ntRe[t] + ii*Global`ntIm[t]] /@ ntTokensOfMonomial[m])]] /@ monos];
  tot = ntIiReduce[Expand[tot], ii];
  {Coefficient[tot, ii, 0], Coefficient[tot, ii, 1]}];

(* Re(plain · Σ_m coef_m Π t_j) = Re(plain)Re(Q) − Im(plain)Im(Q). Splitting it this way is what lets
   `plain` — the dressing coefficient, the denominators, the regulator factors — stay factored: it is
   never multiplied into the expansion, only its two halves are.
   `pureQ` drops every imaginary COEFFICIENT first (the Pure premise Σ Im(c)·tr ≈ 0) but still takes
   the real part of the token PRODUCT, which is exactly what the linear Pure formula got wrong. *)
ntRealOfSummand[s_, pureQ_] := Module[{plain, monos, pr, pi, qr, qi},
  {plain, monos} = ntSplitTokenPart[s];
  {pr, pi} = ntSplitRealImag[plain];
  If[TrueQ[pureQ], pi = 0];
  If[monos === {}, Return[pr]];
  {qr, qi} = ntTokenProductRealImag[monos, pureQ];
  pr*qr - pi*qi];

(* Route by degree, after one shared guard. Only the degree->=2 summands go through the expansion
   above; the linear ones keep their own path (below), which reproduces the projection this file has
   always emitted — same token ordering, same expression shape — so committed kernels stay
   byte-identical. *)
ntProjectIntegrand[integrand_, pureQ_, linear_] := Module[{sums, unsafe, degs, res},
  sums = ntSummandsOf[integrand];
(* GUARD (both projections, before either touches a coefficient). ntPureLinear does not go through
   ntSplitRealImag and has the same blind spot — `(… /. Complex[a_, b_] :> a)` zeroes the I inside
   a + I b just as thoroughly — so this has to sit above the split, not inside it. *)
  unsafe = DeleteDuplicates @ Flatten[
      Function[s, ntIiUnsafeParts[Times @@ Select[ntFactorsOf[s], FreeQ[#, _String] &]]] /@ sums];
  If[unsafe =!= {} && ! TrueQ[$ntComplexRuntimeProjection],
    Message[MakeNTKernel::cplxnest, Length[unsafe], Short[unsafe, 6]]; Abort[]];
  degs = ntTokenDegree /@ sums;
  res =
    If[unsafe =!= {},
(* Runtime projection: EVERY summand takes the exact per-summand path, whatever its token degree.
   The linear fast paths are both unusable here — ntPureLinear's `/. Complex[a_, b_] :> a` deletes
   the I inside a denominator, and ntRePartLinear reaches the same broken extraction through
   ntSplitRealImag. ntRealOfSummand is exact at every degree, including 0 and 1. *)
      ntLog["[complex] ", Length[unsafe], " coefficient(s) carry a Complex below a non-traversable ",
        "head (finite density?); projecting those at generated-code runtime via ntRe/ntIm"];
      Total[ntRealOfSummand[#, pureQ] & /@ sums],
    If[Max[Append[degs, 0]] <= 1,
      linear[integrand],
      Module[{nonlin = Pick[sums, degs, _?(# >= 2 &)]},
        ntLog["[complex] integrand is MULTILINEAR in the trace tokens: ", Length[nonlin], " of ",
          Length[sums], " summand(s), max degree ", Max[degs],
          " (disconnected diagram(s) — exact re/im expansion applied)"];
        linear[Total[Pick[sums, degs, _?(# <= 1 &)]]] + Total[ntRealOfSummand[#, pureQ] & /@ nonlin]]]];
(* No scoped symbol may survive, on EITHER route. This cannot fire once the guard and the routing
   above are right; it is here because the failure mode it guards is a bare identifier in 300 KB of
   C++, diagnosed 20 minutes later by the compiler. (It used to sit after an early Return that the
   all-linear route took, i.e. it never ran on the common case.) See also the $-symbol pattern in
   $ntCppLeakPatterns. *)
  With[{leaked = DeleteDuplicates @ Cases[res, s_Symbol /; StringContainsQ[SymbolName[s], "$"], {0, Infinity}]},
    If[leaked =!= {}, Message[MakeNTKernel::tokleak, leaked]; Abort[]]];
  res];

ntPureIntegrand[integrand_] := ntProjectIntegrand[integrand, True, ntPureLinear];

(* "RePart": the value is real but a trace is itself complex, so only `.real()` of the full complex
   result is correct: Re(Σ c·tr) = Σ[Re(c)·ntRe(tr) − Im(c)·ntIm(tr)]. This is the LINEAR case (each
   term carries one trace token); ntProjectIntegrand routes products of tokens to the multilinear
   split above. The naive `(… /. s:>ntRe[s]+I ntIm[s]) /. Complex[a_,b_]:>a` is WRONG:
   Mathematica keeps `i·X·(ntRe+I ntIm)` as an UNEXPANDED product, so `/. Complex:>a` zeroes the
   leading `i` factor and DROPS the whole term — silently losing the real −X·ntIm(tr) contribution of
   every complex trace. Instead split each token's coefficient into real / imaginary parts via an
   `ii`-substitution (I → real symbol ii), which keeps the dressing coefficients FACTORED (unlike
   ComplexExpand, which un-factors them and defeats COEN's CSE).

   This is done in ONE pass over the summands. The obvious reading of the formula — mint a
   placeholder per token and call Coefficient[integrand, placeholder] once per token — traverses the
   whole integrand once per trace token (three times over, counting the substitution in and the
   token-free remainder out), so it cost `#tokens x expression size`. Bucketing the summands by
   their token instead is flat in the token count and, because Plus is canonically ordered, rebuilds
   the SAME expression: verified `===` identical at 20 / 77 / 154 tokens, so kernels are
   byte-identical.

   HONEST TRADE, measured over 4000 summands (per-token Coefficient -> this):
       5 tokens  0.060 -> 0.105 s      77 tokens  0.138 -> 0.116 s
      20 tokens  0.076 -> 0.105 s     154 tokens  0.235 -> 0.114 s
   i.e. it is ~40% SLOWER below a crossover near 60 tokens and ~2x faster above it. The old shape
   ran a few coarse kernel-internal passes over one big expression; this one runs many fine-grained
   Mathematica-level operations (per summand: ~44 ms in ntLinearSplit, ~56 ms in ntSplitRealImag,
   3 ms in the assembly). Fusing the two helpers would recover ~30 ms of that, at the price of
   duplicating the coefficient split — not taken: 50 ms sits inside a generation that spends minutes
   in NumTrace and net-building, and the reason for the rewrite is not the linear-path speed anyway.

   What it is for: it removes the Unique["tr$"] placeholders from this route entirely — they are
   what leaked into the emitted C++ when the linearity precondition turned out not to hold — and
   leaves ntSplitRealImag as the single place where a coefficient is split, so the ii-substitution
   and its guard (ntIiSafeQ) cannot go out of sync between the two projections. *)

(* One summand -> {its trace token or None, its string-free coefficient}. The caller has already
   established token degree <= 1 for this summand. *)
ntLinearSplit[s_] := Module[{facs = ntFactorsOf[s], tk, rest},
  tk = Cases[facs, _String];
  rest = DeleteCases[facs, _String];
(* The token has to be a BARE factor. One sitting inside a Power or below any other head is invisible
   both to this scan and to the Coefficient extraction it replaces, so such a summand was previously
   dropped outright — a missing term rather than a diagnostic. FreeQ on the leftovers rather than a
   Count over the summand: it short-circuits on the first hit, and this runs once per summand. *)
  If[Length[tk] > 1 || ! FreeQ[rest, _String],
    Message[MakeNTKernel::toknest, s]; Abort[]];
  {If[tk === {}, None, First[tk]], Times @@ rest}];

ntRePartLinear[integrand_] := Module[{parts, groups, toks, tokPart, constPart},
  parts = Function[s, Module[{tk, plain}, {tk, plain} = ntLinearSplit[s]; {tk, ntSplitRealImag[plain]}]] /@
      ntSummandsOf[integrand];
  groups = GroupBy[parts, First -> Last];
(* Sort matches the Union[Cases[…]] ordering the per-token formulation emitted in. It is presentation
   only — Plus is Orderless, so Total canonicalises regardless, and reversing this line was verified
   to leave the result `===` unchanged. Kept so a printed intermediate reads in token order. *)
  toks = Sort[DeleteCases[Keys[groups], None]];
  tokPart = Total[Function[t, With[{v = groups[t]},
       Total[v[[All, 1]]]*Global`ntRe[t] - Total[v[[All, 2]]]*Global`ntIm[t]]] /@ toks];
  constPart = If[KeyExistsQ[groups, None], Total[groups[None][[All, 1]]], 0];
  tokPart + constPart];

ntRePartIntegrand[integrand_] := ntProjectIntegrand[integrand, False, ntRePartLinear];
