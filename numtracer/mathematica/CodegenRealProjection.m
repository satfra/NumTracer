(* Code generation: the two REAL projections of a complex integrand (Pure and RePart), including the
   exact multilinear re/im split and the guard against non-polynomial Complex coefficients.
   Loaded by NumTracer.m via ntLoadPart, in the NumTracer`Private` context. *)

(* ---- the two REAL projections of a complex integrand ---------------------------------------
   Both are emitted for every complex flow; the imaginary-part probe picks between them (and the
   untouched complex form) with a preprocessor `#if` — see ntProbeSource (CodegenProbe.m).

   "Pure": the projection Complex -> Re is exact, i.e. Σ Im(c)·tr ≈ 0. Wrap the trace tokens in ntRe
   FIRST so the kernel is provably double-typed even when a trace function is complex-typed, then drop
   the imaginary coefficients outright. *)
ntPureLinear[integrand_] := (integrand /. s_String :> Global`ntRe[s]) /. Complex[a_, b_] :> a;

(* ---- token degree: the multilinear generalisation -------------------------------------------
   The integrand is NOT always linear in the trace tokens. A disconnected diagram (factorIdsOf =!= None)
   contributes `Π traceRef[factor groups] · traceRef[anchor]`, a PRODUCT of tokens (see
   ntAssembleIntegrand). The linear projections are wrong on such summands: Coefficient leaves the
   other tokens inside, a term is counted once per token, t^2 is dropped, and ntPureLinear's
   per-token ntRe loses the −Im(A)·Im(B) leg of Re(c·A·B) (compiles, wrong number).

   Exact for any degree, per SUMMAND: keep the string-free part FACTORED (for COEN's CSE; this is why
   ComplexExpand is not used), expand only the token-bearing part (a product of short sums), and take
   the real part of each monomial with the `ii`-substitution. Cost is 2^n products at degree n. *)

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
ntReduceImagUnitPowers[e_, ii_] := e //. Power[ii, n_Integer /; n >= 2] :> (-1)^Quotient[n, 2] * ii^Mod[n, 2];

(* ---- when is the ii-substitution faithful? --------------------------------------------------
   Substituting I -> ii (a real symbol) and reading off Coefficient[·, ii, 0|1] is EXACT only where
   the expression is a POLYNOMIAL in ii. Otherwise Mathematica silently series-expands and returns the
   leading term:
     "tr0"/(a + I b)^n ->  ntRe["tr0"]/a^n      imaginary part of the denominator deleted (at finite
                                                 mu, muq vanishes from every quark propagator)
     Sqrt[a + I b]     ->  Sqrt[a + b ii$7182]  the Module-local stand-in leaks into the C++
   Such shapes are also what made the projection slow (tens of seconds at a few hundred summands).

   The test is purely STRUCTURAL (no Expand, so it is free on large factored coefficients). It admits
   a Complex, products and sums of admissible parts, and a non-negative integer power of one;
   everything else carrying a Complex is refused. *)
ntImagSplitSafeQ[e_] :=
  FreeQ[e, Complex] ||
  Switch[Head[e],
    Complex, True,
    Times | Plus, AllTrue[List @@ e, ntImagSplitSafeQ],
    Power, IntegerQ[e[[2]]] && NonNegative[e[[2]]] && ntImagSplitSafeQ[e[[1]]],
    _, False];

(* The MAXIMAL offending nodes, for the message: descend through the Times/Plus skeleton and stop at
   the first inadmissible node, so the message shows the offending denominator, not the whole summand. *)
ntImagSplitUnsafeParts[e_] :=
  If[ntImagSplitSafeQ[e], {},
    Switch[Head[e],
      Times | Plus, Flatten[ntImagSplitUnsafeParts /@ (List @@ e)],
      _, {e}]];

(* {Re[e], Im[e]} for a string-free coefficient whose only imaginary content is explicit `Complex`
   numbers (no symbol here is imaginary). The common shapes — real, or one Complex scaling a factored
   real expression — are taken WITHOUT expanding, which keeps dressing coefficients factored for
   COEN's CSE. Anything else falls back to the ii-substitution. *)
ntSplitRealImag[e_] := Module[{fs, cs, rest, c, ii, t},
  If[FreeQ[e, Complex], Return[{e, 0}]];
  fs = ntFactorsOf[e];
  cs = Cases[fs, _Complex];
  rest = DeleteCases[fs, _Complex];
  If[cs =!= {} && FreeQ[rest, Complex],
    c = Times @@ cs;
    Return[{Re[c] * Times @@ rest, Im[c] * Times @@ rest}]];
  (* Not a polynomial in the stand-in (see ntImagSplitSafeQ): under "ComplexRuntimeProjection" hand the
     untouched, still-factored expression to the generated code to split at runtime. Otherwise
     ntProjectIntegrand has already refused it, so the extraction below only sees safe input. *)
  If[! ntImagSplitSafeQ[e] && TrueQ[$ntComplexRuntimeProjection],
    Return[{Global`ntRe[e], Global`ntIm[e]}]];
  t = ntReduceImagUnitPowers[Expand[e /. Complex[ar_, ai_] :> ar + ii*ai], ii];
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
  tot = ntReduceImagUnitPowers[Expand[tot], ii];
  {Coefficient[tot, ii, 0], Coefficient[tot, ii, 1]}];

(* Re(plain · Σ_m coef_m Π t_j) = Re(plain)Re(Q) − Im(plain)Im(Q). Splitting it this way keeps
   `plain` (dressing, denominators, regulators) factored: only its two halves enter, never the expansion.
   `pureQ` drops every imaginary COEFFICIENT (the Pure premise Σ Im(c)·tr ≈ 0) but still takes the
   real part of the token PRODUCT. *)
ntRealOfSummand[s_, pureQ_] := Module[{plain, monos, pr, pi, qr, qi},
  {plain, monos} = ntSplitTokenPart[s];
  {pr, pi} = ntSplitRealImag[plain];
  If[TrueQ[pureQ], pi = 0];
  If[monos === {}, Return[pr]];
  {qr, qi} = ntTokenProductRealImag[monos, pureQ];
  pr*qr - pi*qi];

(* NB: no backquoted code fragments in message strings: a backquoted word is a StringForm SLOT, so
   quoting an identifier that way makes the message itself fail to format (StringForm::sfr). *)

MakeNTKernel::cplxnest = "ntProjectIntegrand: a Complex sits below a head the real/imaginary split cannot traverse (typically a denominator such as the finite-density l0 + I muq). The split is exact only for coefficients polynomial in I; here it would silently drop the imaginary part of the denominator. Pass \"ComplexRuntimeProjection\" -> True to project at runtime instead. `1` offending subexpression(s):\n`2`";

MakeNTKernel::tokleak = "ntProjectIntegrand: a scoped symbol survived the real/imaginary projection and would be emitted as a bare C++ identifier (which compiles). Either the stand-in for I survived a coefficient that is not polynomial in it, or a trace-token placeholder was not substituted back: an integrand shape the projection routing does not cover. Offending symbol(s):\n`1`";

(* Route by degree, after one shared guard. Only degree >= 2 summands go through the expansion above;
   the linear ones keep the linear path, whose expression shape keeps committed kernels byte-identical. *)
ntProjectIntegrand[integrand_, pureQ_, linear_] := Module[{sums, unsafe, degs, res},
  sums = ntSummandsOf[integrand];
  (* GUARD for both projections, before either touches a coefficient. It sits here rather than in
     ntSplitRealImag because ntPureLinear bypasses that and has the same blind spot: its
     `/. Complex[a_, b_] :> a` zeroes the I inside a + I b. *)
  unsafe = DeleteDuplicates @ Flatten[
      Function[s, ntImagSplitUnsafeParts[Times @@ Select[ntFactorsOf[s], FreeQ[#, _String] &]]] /@ sums];
  If[unsafe =!= {} && ! TrueQ[$ntComplexRuntimeProjection],
    Message[MakeNTKernel::cplxnest, Length[unsafe], Short[unsafe, 6]]; Abort[]];
  degs = ntTokenDegree /@ sums;
  res =
    If[unsafe =!= {},
      (* Runtime projection: EVERY summand takes the exact per-summand path, whatever its degree. Both
         linear paths would mangle the denominators (ntPureLinear directly, ntRePartLinear via the
         ii-extraction); ntRealOfSummand is exact at every degree, including 0 and 1. *)
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
  (* No scoped symbol may survive, on either route. Should be unreachable given the guard and routing
     above; it exists because the failure is otherwise a bare identifier in the C++, found only by the
     compiler. See also the $-symbol pattern in $ntCppLeakPatterns. *)
  With[{leaked = DeleteDuplicates @ Cases[res, s_Symbol /; StringContainsQ[SymbolName[s], "$"], {0, Infinity}]},
    If[leaked =!= {}, Message[MakeNTKernel::tokleak, leaked]; Abort[]]];
  res];

ntPureIntegrand[integrand_] := ntProjectIntegrand[integrand, True, ntPureLinear];

(* "RePart": the value is real but a trace is itself complex, so only `.real()` of the full complex
   result is correct: Re(Σ c·tr) = Σ[Re(c)·ntRe(tr) − Im(c)·ntIm(tr)]. This is the LINEAR case; products
   of tokens are routed to the multilinear split by ntProjectIntegrand.
   TRAP: the naive `(… /. s:>ntRe[s]+I ntIm[s]) /. Complex[a_,b_]:>a` is WRONG. Mathematica keeps
   `i·X·(ntRe+I ntIm)` as an unexpanded product, so `/. Complex:>a` zeroes the leading `i` and drops
   the real −X·ntIm(tr) contribution of every complex trace. Hence the explicit coefficient split.

   One pass: summands are bucketed by their token instead of calling Coefficient once per token (cost
   #tokens × expression size). Plus is canonically ordered, so the result is `===` the per-token form.
   ntSplitRealImag stays the single place a coefficient is split, so the ii-substitution and its guard
   (ntImagSplitSafeQ) cannot diverge between the two projections. *)

(* One summand -> {its trace token or None, its string-free coefficient}. The caller has already
   established token degree <= 1 for this summand. *)
ntLinearSplit[s_] := Module[{facs = ntFactorsOf[s], tk, rest},
  tk = Cases[facs, _String];
  rest = DeleteCases[facs, _String];
  (* The token has to be a BARE factor: one inside a Power or another head is invisible to this scan,
     and the summand would be dropped silently. Abort instead. *)
  If[Length[tk] > 1 || ! FreeQ[rest, _String],
    Message[MakeNTKernel::toknest, s]; Abort[]];
  {If[tk === {}, None, First[tk]], Times @@ rest}];

MakeNTKernel::toknest = "ntRePartLinear: a trace token in this summand is not a BARE factor of it — it sits inside a Power, or below some other head. The token-degree routing classified the summand as linear, but neither the factor-level scan here nor the Coefficient extraction it replaced can see such a token, so the whole summand would be dropped from the projected integrand: a missing term, with nothing downstream to notice. Offending summand:\n`1`";

ntRePartLinear[integrand_] := Module[{parts, groups, toks, tokPart, constPart},
  parts = Function[s, Module[{tk, plain}, {tk, plain} = ntLinearSplit[s]; {tk, ntSplitRealImag[plain]}]] /@
      ntSummandsOf[integrand];
  groups = GroupBy[parts, First -> Last];
  (* Sort is presentation only (Plus is Orderless); it makes a printed intermediate read in token order. *)
  toks = Sort[DeleteCases[Keys[groups], None]];
  tokPart = Total[Function[t, With[{v = groups[t]},
       Total[v[[All, 1]]]*Global`ntRe[t] - Total[v[[All, 2]]]*Global`ntIm[t]]] /@ toks];
  constPart = If[KeyExistsQ[groups, None], Total[groups[None][[All, 1]]], 0];
  tokPart + constPart];

ntRePartIntegrand[integrand_] := ntProjectIntegrand[integrand, False, ntRePartLinear];
