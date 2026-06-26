(* ::Package:: *)

(* Code generation. NumTracer owns ONLY the tensor part — turning each contraction
   component into an `et` builder chain whose pruned scalar is read out by
   expr::eval_real. Everything scalar (the dressing/regulator coefficients, CSE,
   powr<n>, the function/class/header boilerplate, clang-format, write-if-changed)
   is delegated to FunKit's mature COEN emitter (CppForm / MakeCppFunction /
   MakeCppClass / MakeCppHeader / WriteCodeToFile), which produces the flat
   straight-line kernel form.

   The seam: each component's scalar result is bound to a C++ identifier in the
   kernel-body PREAMBLE we emit; that identifier appears as a *string placeholder*
   in the Mathematica integrand `Σ coeff_i × trace_i`, which FunKit's CppForm emits
   verbatim (CExpression[a_String] := a) while CSE-ing the coefficients around it. *)

(* The symbolic ET expression-tensor backend (`builder`/`compileT`/`compileDiagram`, emitting
   `numtracer::et::` calls) has been REMOVED. Generation goes through the numeric matrix-product
   backend (`builderInv`/`compileTInv` + the generator below); the Dense `numtracer::dense::` oracle
   path remains for cross-checking. *)

(* ============================================================================ *)
(* ==== DENSE path: the brute-force numeric backend. Same tensor structure as  === *)
(* ==== compileT, but emitting `numtracer::dense::` DTensor builder calls that  === *)
(* ==== contract NUMBERS (all index combinations, no structural-zero pruning).  === *)
(* ==== The momentum builders read the same `double renv[]` the et path fills;  === *)
(* ==== a component's scalar is read directly via `.scalar_value().re` (constexpr *)
(* ==== for a colour/gamma constant, else a runtime const). This is the naive    *)
(* ==== dense-matrix baseline benchmarked against the pruned et / reduced inv.    *)
(* ============================================================================ *)

builderDense[ntMetric[mu_, nu_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::metric<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ">()";
builderDense[ntVec[q_, mu_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::vector<" <> ToString[ids[mu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ">(renv)";
builderDense[ntTransProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::transverse_projector<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <>
    ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ">(renv)";
builderDense[ntLongProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::longitudinal_projector<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <>
    ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ">(renv)";
builderDense[ntMagneticProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::magnetic_projector<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <>
    ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["InvS"]] <> ">(renv)";
builderDense[ntElectricProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::electric_projector<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <>
    ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <>
    ", " <> ToString[env[q]["InvS"]] <> ">(renv)";
builderDense[ntSUNf[n_, a_, b_, c_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::f<" <> ToString[n] <> ", " <> ToString[ids[a]] <> ", " <> ToString[ids[b]] <> ", " <> ToString[ids[c]] <> ">()";
builderDense[ntSUNDeltaAdj[n_, a_, b_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::delta_adj<" <> ToString[n] <> ", " <> ToString[ids[a]] <> ", " <> ToString[ids[b]] <> ">()";
builderDense[ntGamma[mu_, din_, dout_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::gamma_axis<" <> ToString[ids[mu]] <> ", " <> ToString[ids[din]] <> ", " <> ToString[ids[dout]] <> ">()";
builderDense[ntGamma5[din_, dout_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::gamma5<" <> ToString[ids[din]] <> ", " <> ToString[ids[dout]] <> ">()";
builderDense[ntDeltaDirac[din_, dout_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::identity<" <> ToString[ids[din]] <> ", " <> ToString[ids[dout]] <> ">()";
builderDense[ntSUNT[n_, a_, i_, j_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::T<" <> ToString[n] <> ", " <> ToString[ids[a]] <> ", " <> ToString[ids[i]] <> ", " <> ToString[ids[j]] <> ">()";
builderDense[ntSUNDeltaFund[n_, i_, j_], ids_, env_, mask_, nc_] :=
  "numtracer::dense::delta_fund<" <> ToString[n] <> ", " <> ToString[ids[i]] <> ", " <> ToString[ids[j]] <> ">()";

scaleStrDense[str_, 1] := str;
scaleStrDense[str_, 1.] := str;
scaleStrDense[str_, s_] := "numtracer::dense::scale<numtracer::Cx{" <> ToString[CForm[N[s, 17]]] <> ", 0.0}>(" <> str <> ")";
wrapContractDense[{one_}] := one;
wrapContractDense[many_] := "numtracer::dense::contract_all(" <> StringRiffle[many, ", "] <> ")";

compileTDense[e_, ids_, env_, mask_, nc_] := Which[
  tensorQ[e],
    {builderDense[e, ids, env, mask, nc], 1},
  Head[e] === Times,
    Module[{parts = List @@ e, sc, tn, cs},
      sc = Select[parts, scalarQ]; tn = Select[parts, ! scalarQ[#] &];
      cs = compileTDense[#, ids, env, mask, nc] & /@ orderFactors[tn];
      {wrapContractDense[cs[[All, 1]]], (Times @@ sc) (Times @@ cs[[All, 2]])}],
  Head[e] === Plus,
    Module[{cs = compileTDense[#, ids, env, mask, nc] & /@ (List @@ e)},
      If[! AllTrue[cs[[All, 2]], NumericQ], Message[MakeNTKernel::eagernn, e]; Abort[]];
      {"numtracer::dense::add_all(" <> StringRiffle[MapThread[scaleStrDense, {cs[[All, 1]], cs[[All, 2]]}], ", "] <> ")", 1}],
  True, {"", e}];

(* One diagram, dense backend. Same shape as compileDiagram, but the scalar is read
   straight off the numeric DTensor — no escalar_t/eval_real. When `cplx` (the kernel has
   imaginary projector/propagator coefficients), the momentum-dependent trace tokens keep
   their FULL complex value (`.scalar_cplx()`); colour/gamma constants stay real doubles, and
   the kernel takes the real part of the assembled integrand. *)
compileDiagramDense[diag_, tag_, frame_, env_, mask_, nc_, cplx_] := Module[{ids = diag["Ids"], coeff, decls = {}, toks = {}},
  coeff = diag["Coeff"] /. {
    ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame],
    ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]],
    ntVec[q_, i_Integer] :> resolveComponents[q, frame][[i + 1]]};
  MapIndexed[Function[{comp, ci},
    Module[{str, scal, ct, trVar},
      {str, scal} = compileTDense[Times @@ comp["Factors"], ids, env, mask, nc];
      coeff *= scal;
      ct = tag <> "_" <> ToString[ci[[1]] - 1];  trVar = "_tr" <> ct;
      decls = Join[decls, {
        Which[
          (* complex kernel: keep EVERY component complex. The colour factor of a non-abelian
             vertex diagram (f^{abc} T^b T^c = (iN/2) T^a) is IMAGINARY — taking its real part
             drops the diagram entirely. Its i pairs with the projector's i to give a real flow. *)
          cplx && comp["Constant"], "const std::complex<double> " <> trVar <> " = " <> str <> ".scalar_cplx();",
          comp["Constant"],         "constexpr double " <> trVar <> " = " <> str <> ".scalar_value().re;",
          cplx,                     "const std::complex<double> " <> trVar <> " = " <> str <> ".scalar_cplx();",
          True,                     "const double " <> trVar <> " = " <> str <> ".scalar_value().re;"]}];
      toks = Append[toks, trVar]]], diag["Components"]];
  {decls, coeff, toks}
];

(* ============================================================================ *)
(* ==== INVARIANT-BASIS path: contract in scalar-product symbols, not frame  === *)
(* ==== components. The Lorentz trace of each diagram is emitted as an `et::inv` *)
(* ==== network into a build-time GENERATOR (gen_<name>_inv.cpp); the generator *)
(* ==== runs reduce->cancel->lower at codegen time (heap reclaimed, ~FORM op    *)
(* ==== count) and prints a committed straight-line `<name>_inv_kernels.hh`.    *)
(* ==== The kernel just fills the few fundamental symbols and calls trN(f).     *)
(* ==== This closes the A3 runtime gap that the frame-component basis cannot.   *)
(* ============================================================================ *)

(* Verbose-diagnostics gate. The profiling / CSE / probe / timing traces below ([prof], [cse],
   [probe], [diagpoly], [time]) are emitted through ntLog and stay SILENT unless this flag is set —
   so a normal generation run is quiet. Genuine "[NumTracer] ERROR" aborts and the "wrote:"/
   "unchanged:" file messages are always printed (plain Print). To see the diagnostics, set
   NumTracer`Private`$NumTracerVerbose = True before generating. ntLog evaluates its arguments only
   when verbose, so it must NEVER wrap a side-effecting computation — keep such work (e.g. an
   AbsoluteTiming around an assignment) in a With-binding and pass only the timing into ntLog. *)
$NumTracerVerbose = False;
ntLog[args___] := If[TrueQ[$NumTracerVerbose], Print[args]];

(* a C++ double literal at full precision (exact rationals -> doubles). *)
cppNum[x_] := ToString[CForm[N[x, 17]]];

(* CppForm, flattened to a single line — the fill formulas are emitted INSIDE C++ string
   literals, so CppForm's line-wrapping (newlines) must be collapsed or they break the string. *)
cppFlat[e_] := StringReplace[FunKit`CppForm[e], {"\n" -> " ", "\r" -> " ", "\t" -> " "}];

(* ---- GlobalCollect: dressing-coefficient decomposition (the Route-B front-end) ----------------
   A diagram's dressing/kinematic coefficient is decomposed into a SUM of monomials, each
   {numericCoeff, {atomExpr...}}. An ATOM is a maximal non-numeric multiplicative factor that the
   kernel evaluates ONCE (a propagator dressing `ZA3[...]`, a regulator `RB[...]`, a composite
   denominator `Power[D,-1]`, a kinematic factor `cos1`/`Sqrt[...]`). Power[b,n] with a positive
   integer n expands to n copies of atom b (so it lowers as b^n); any other factor — including a
   negative/fractional power or a function call — is one atom. The numeric part folds into the
   collector coefficient; the atoms become inert DRESSING SYMBOLS (a `dr` tag) so identical
   (scalar-product × dressing) monomials from different diagrams MERGE — the cross-diagram collection
   FORM does. Returns {monomials, ...} where each monomial is {num (machine real), {atomExpr...}}. *)
dressMonomials[expr_] := Module[{terms},
  terms = Expand[expr];
  terms = If[Head[terms] === Plus, List @@ terms, {terms}];
  Function[term, Module[{num = 1, atoms = {}, factors},
    factors = If[Head[term] === Times, List @@ term, {term}];
    Do[Which[
        NumberQ[f], num *= f,
        MatchQ[f, Power[_, _Integer?Positive]], atoms = Join[atoms, ConstantArray[First[f], Last[f]]],
        True, AppendTo[atoms, f]],
      {f, factors}];
    {N[num, 17], atoms}]] /@ terms];

(* ---- inv dialect: same tensor structure as `compileT`, but emitting the
        et::inv builder calls the generator's wrapper templates expose
        (tproj / lmetric / lvec / sc / contract / add). *)
builderInv[ntMetric[mu_, nu_], ids_, env_, mask_, nc_] :=
  "lmetric<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ">()";
builderInv[ntVec[q_, mu_], ids_, env_, mask_, nc_] :=
  "lvec<" <> ToString[ids[mu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ">()";
builderInv[ntTransProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "tproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <>
    ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ">()";
builderInv[ntLongProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "lproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <>
    ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ">()";
builderInv[ntMagneticProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "mproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <>
    ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["InvS"]] <> ">()";
builderInv[ntElectricProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "eproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <>
    ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ", " <> ToString[env[q]["InvS"]] <> ">()";
builderInv[ntEpsilon[a_, b_, c_, d_], ids_, env_, mask_, nc_] :=
  "leps<" <> ToString[ids[a]] <> ", " <> ToString[ids[b]] <> ", " <> ToString[ids[c]] <> ", " <>
    ToString[ids[d]] <> ">()";

scaleStrInv[str_, 1] := str;
scaleStrInv[str_, 1.] := str;
scaleStrInv[str_, s_] := "sc<numtracer::expr::Lit<numtracer::Cx{" <> cppNum[s] <> ", 0.0}>>(" <> str <> ")";
wrapContractInv[{one_}] := one;
wrapContractInv[many_] := "contract(" <> StringRiffle[many, ", "] <> ")";

(* MEMOIZED: the projector/Lorentz net builder is called with ~95% repeated inputs (the same
   transverse-projector structures recur across every diagram/branch — measured 166296 calls but only
   7264 distinct on ZA3 1/4/7). Cache by a hash of the args; $ctCache is cleared per generation in
   mkGenerateKernel. Output-preserving (a pure function of its inputs); the recursion is memoised too. *)
compileTInv[e_, ids_, env_, mask_, nc_] := With[{h = Hash[{e, ids, env, mask, nc}]},
  Lookup[$ctCache, h, $ctCache[h] = compileTInvBody[e, ids, env, mask, nc]]];
compileTInvBody[e_, ids_, env_, mask_, nc_] := Which[
  tensorQ[e],
    {builderInv[e, ids, env, mask, nc], 1},
  Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 1 && ! scalarQ[e],
    (* A TENSOR raised to an integer power = that many copies sharing the SAME index labels, i.e. a
       closed self-contraction (e.g. ntMetric[v1,v2]^2 = g_{v1 v2} g_{v1 v2} = D). Expand into n
       contracted copies so the numeric index-elimination folds it to a number. Without this it falls
       through to the scalar branch below and is CForm'd into undeclared C++ (the ZAAqbq metric leak). *)
    Module[{cs = Table[compileTInv[e[[1]], ids, env, mask, nc], {e[[2]]}]},
      {wrapContractInv[cs[[All, 1]]], Times @@ cs[[All, 2]]}],
  Head[e] === Times,
    Module[{parts = List @@ e, sc, tn, cs},
      sc = Select[parts, scalarQ]; tn = Select[parts, ! scalarQ[#] &];
      cs = compileTInv[#, ids, env, mask, nc] & /@ orderFactors[tn];
      {wrapContractInv[cs[[All, 1]]], (Times @@ sc) (Times @@ cs[[All, 2]])}],
  Head[e] === Plus,
    Module[{cs = compileTInv[#, ids, env, mask, nc] & /@ (List @@ e)},
      If[! AllTrue[cs[[All, 2]], NumericQ], Message[MakeNTKernel::eagernn, e]; Abort[]];
      {"add(" <> StringRiffle[MapThread[scaleStrInv, {cs[[All, 1]], cs[[All, 2]]}], ", "] <> ")", 1}],
  True, {"", e}];

(* ---- colour/Lorentz sector split for a colour-ENTANGLED Lorentz/Dirac component -----------
   The full quark-gluon vertex basis (AqbqDirect147 structures 4/7) produces a quark-loop
   component that is a SUM whose terms pair different colour orderings (T^{c1}T^{c2} vs
   T^{c2}T^{c1}) with different Dirac traces — colour and Dirac do NOT factor globally, only per
   term. expandDiracComponent then leaves the colour tensors inside the (would-be Lorentz) result,
   where the Lorentz-only builderInv cannot lower them. The fix: group the Dirac-traced expression
   by its colour-factor product; each group has ONE colour structure (folded numerically as a
   SUNNet — exactly like a constant adjoint/fundamental component) times a pure-Lorentz
   polynomial (the inv net). The diagram's trace is Sum_group colv_group * lorentzPoly_group,
   emitted as several (sunNet, lorentzNet) entries that the per-diagram combination already sums. *)
ctHeadsInv = {_ntSUNf, _ntSUNDeltaAdj, _ntSUNT, _ntSUNDeltaFund, _ntSUNDiagFund, _ntSUNDiagAdj};
colourEntangledQ[e_] := ! FreeQ[e, Alternatives @@ ctHeadsInv];
colourFactorProd[term_] := Times @@ Cases[If[Head[term] === Times, List @@ term, {term}], Alternatives @@ ctHeadsInv];
mergeColNet["SUNNet{}", b_] := b;
mergeColNet[a_, "SUNNet{}"] := a;
mergeColNet[a_, b_] := "SUNNet{" <> StringDrop[StringDrop[a, 7], -1] <> "," <> StringDrop[StringDrop[b, 7], -1] <> "}";

(* ---- per-component diagonal-dressing registry (ntSUNDiag{Fund,Adj}) ------------------------
   A diag-dressed group δ carries a runtime per-component dressing `name` (a kernel std::array)
   evaluated at a kinematic `scale`. colFacG registers each distinct (name, scale) under a small
   integer dressing id `dr` (baked into the emitted SUN::diag{Fund,Adj}(...,dr) factor). The C++
   seam then folds the net to a SUNPoly over these ids, and the integrand multiplies in the runtime
   sum Σ_t coeff_t Π name[comp](scale) via the ntDiagDress(arr,comp,scale) -> arr[comp](scale)
   helper. Reset per generation. *)
$diagDrTable = <||>; $diagDrByKey = <||>; $diagDrCounter = 0;
resetDiagDr[] := ($diagDrTable = <||>; $diagDrByKey = <||>; $diagDrCounter = 0;);
diagDrId[name_, scale_, dim_, kind_] := Module[{key = {name, scale}},
  If[! KeyExistsQ[$diagDrByKey, key],
    $diagDrByKey[key] = $diagDrCounter;
    $diagDrTable[$diagDrCounter] = <|"Name" -> name, "Scale" -> scale, "Dim" -> dim, "Kind" -> kind|>;
    $diagDrCounter++];
  $diagDrByKey[key]];

(* ---- scalar-dressing registry (symbolic dressing collection: ntDressedNum slots) -------------
   Each DISTINCT dressing atom — a maximal non-numeric multiplicative factor in a dressed numerator's
   per-structure coefficient (a propagator dressing `Zq[...]`, `hSigL[...]`, a regulator `RB[...]`, a
   composite denominator `Power[D,-1]`, a kinematic `cos1`) — is interned under a small integer id baked
   into the emitted DSlotOpt and into the generator's fm.dress table. So identical atoms across slots /
   diagrams share ONE `f[]` slot and the runtime evaluates each dressing call ONCE (the cross-diagram
   collection FORM does). The atom expression is already frame-resolved (ntSP/ntVec components
   substituted) so its C++ fill is `cppFlat[atom]`. Reset per generation alongside resetDiagDr. *)
$drTable = <||>; $drByKey = <||>; $drCounter = 0;
resetDr[] := ($drTable = <||>; $drByKey = <||>; $drCounter = 0;);
drAtomId[atom_] := Module[{key = atom},
  If[! KeyExistsQ[$drByKey, key],
    $drByKey[key] = $drCounter; $drTable[$drCounter] = atom; $drCounter++];
  $drByKey[key]];
(* Decompose a (frame-resolved) dressed-structure coefficient into {Cx numeric, {drAtomId…}}: numbers
   fold into the complex coefficient; a positive-integer power b^n expands to n atom copies; every other
   non-numeric factor is one atom (interned via drAtomId). Mirrors dressMonomials but interns atoms. *)
drDecompose[coeff_] := Module[{factors = If[Head[coeff] === Times, List @@ coeff, {coeff}], num = 1, ids = {}},
  Do[Which[
      NumberQ[f], num *= f,
      MatchQ[f, Power[_, _Integer?Positive]], ids = Join[ids, ConstantArray[drAtomId[First[f]], Last[f]]],
      True, AppendTo[ids, drAtomId[f]]], {f, factors}];
  {N[num, 17], Sort[ids]}];

(* Chunk a Lorentz polynomial into several inv nets of <= $ntInvChunk top-level terms each, so no
   single generated net-builder function becomes a giant nested add() that blows up the g++ -O0
   compile (a quark box with the full quark-gluon vertex basis can be tens of thousands of nodes in
   ONE net). Returns a list of {lorentzNetString, scalar}; the per-diagram combination sums the
   chunks (same colour + coeff) into one trace, so chunking is transparent. Zero -> dropped; a pure
   number -> a constant net (konst). *)
$ntInvChunk = 150;
(* max concurrent g++ processes in the generator compile (peak RAM ~ jobs x per-unit). Scales with
   the machine: NT_GEN_JOBS env override, else min(16, cores-2) — the -O0 net-builder TUs are light
   on RAM, and dense flows (e.g. the σ^μν struct-7 projection) emit many large units, so the default
   4 left most cores idle (a 51-unit struct-7 compile ran ~1 h at -P 4). Capped at 16 for RAM. *)
$ntCompileJobs = With[{e = Environment["NT_GEN_JOBS"]},
  Which[
    StringQ[e] && IntegerQ[Quiet@ToExpression[e]] && ToExpression[e] > 0, ToExpression[e],
    IntegerQ[$ProcessorCount], Max[4, Min[16, $ProcessorCount - 2]],
    True, 4]];

(* C++ compiler for the build-time generator. The emitted generator is ordinary numeric C++, so
   either g++ or clang++ compiles it with the same flags (-std=c++20 -ftemplate-depth=4000
   -O2/-O0 -pthread -I -c -o). Resolution order: the NT_GEN_CXX env override (verbatim), else an
   installed compiler detected via CCompilerDriver (GCC -> g++, Clang -> clang++; the binary is
   taken from the driver's CompilerInstallation dir when it exists there), else "g++" on PATH. *)
resolveGenCxx[] := Module[{env = Environment["NT_GEN_CXX"], comps, pick, name, inst, exe, full},
  If[StringQ[env] && StringTrim[env] =!= "", Return[StringTrim[env]]];
  Quiet@Needs["CCompilerDriver`"];
  comps = Quiet@Check[CCompilers[], {}];
  If[! ListQ[comps], comps = {}];
  comps = Select[comps, AssociationQ[#] && StringQ[Lookup[#, "CompilerInstallation"]] &];
  If[comps === {}, Return["g++"]];
  pick = SelectFirst[comps, Lookup[#, "Compiler"] === Quiet@DefaultCCompiler[] &, First[comps]];
  name = ToString@Lookup[pick, "Name", ""];
  inst = Lookup[pick, "CompilerInstallation"];
  exe = If[StringContainsQ[name, "Clang", IgnoreCase -> True], "clang++", "g++"];
  full = FileNameJoin[{inst, exe}];
  If[FileExistsQ[full], full, exe]];
chunkLorInv[lorExpr_, ids_, env_, mask_, nc_] := Which[
  lorExpr === 0, {},
  scalarQ[lorExpr], {{"konst(" <> cppNum[lorExpr] <> ")", 1}},
  True, Module[{terms = If[Head[lorExpr] === Plus, List @@ lorExpr, {lorExpr}]},
    (compileTInv[Total[#], ids, env, mask, nc] & /@ Partition[terms, UpTo[$ntInvChunk]])]];

(* Build inv nets for keepLor * inner, where keepLor is the big COMMON Lorentz polynomial (the
   projector x kinematics, shared by every branch of a colour group) and inner is the per-group
   Dirac-trace sum. keepLor's terms are chunked (<= $ntInvChunk each) and inner is multiplied onto each
   chunk, so the projector net is built ONCE PER CHUNK instead of once per branch. Building the
   projector net per branch (the old `chunkLorInv` over a fully-multiplied-out lorSum) was the
   string-build hot spot — a single 3-/4-point quark component could sit minutes inside compileTInv. *)
chunkLorProd[keepLor_, inner_, ids_, env_, mask_, nc_] := Which[
  inner === 0 || inner === 0., {},
  keepLor === 1, chunkLorInv[inner, ids, env, mask, nc],
  Head[keepLor] === Plus,
    (compileTInv[Total[#] inner, ids, env, mask, nc] & /@ Partition[List @@ keepLor, UpTo[$ntInvChunk]]),
  True, {compileTInv[keepLor inner, ids, env, mask, nc]}];

(* {colourNet, lorentzNet, scalar} entries per colour-structure group of a colour-entangled component.
   DISTRIBUTE only the entangled Pluses (the vertex colour-ordering / commutator sub-sums, small) — NOT
   the big pure-Lorentz angular polynomial, which stays factored in each branch. Each branch is split
   into its colour-factor product (folded numerically → SUNNet) and its colour-free rest; the rest
   is emitted as a `dirac_value` net (gammas traced in C++) when it carries a gamma chain, else as a
   chunked Lorentz net. Branches are grouped by colour structure and summed. The Dirac trace is thus
   contracted in the C++ generator, not expanded symbolically here. *)
(* ---- struct-7 σ^{μν} folding: keep the bare γ-commutator as ONE token ------------------------
   The quark-gluon-vertex struct-7 tensor σ^{μν}=(i/2)[γ^μ,γ^ν] arrives from FunKit as a bare 2-term
   antisymmetric γ-pair `Plus`  s·(γ(X)γ(Y) − γ(Y)γ(X))  (FunKit has no σ primitive). Left alone,
   splitColourGroupsInv's `Expand` distributes it into TWO full Dirac traces. `foldDiracSigma`
   collapses that Plus into a single `ntSigma[legA, legB, din, dout]` token (each leg a slashed
   momentum {"slash",mom} or a free gluon id {"free",mu}), so the commutator is traced ONCE (the C++
   engine folds [A,B] as a block-diagonal 2×2 factor — `dcomm*` in et/inv/dirac.hpp). The i/2 and any
   sign live in the SCALAR (the Plus is already a bare bracket). It is a pure OPTIMIZATION: when the
   Plus is not an UNAMBIGUOUS commutator the recognizer returns $Failed and it distributes as before. *)

(* the (type,value) leg of a gamma `g` within a term's factor list `tf`. The gamma's Lorentz label μ
   is carried by the factor(s) holding ntVec[_,μ): either a single ntVec[mom,μ] OR a nested momentum
   sum like (ntVec[l1,μ]−ntVec[p2,μ]) — the LOOP σ legs are momentum linear combinations (l1−p2, …),
   not single vecs. The leg is returned as a SORTED list of {coeff, q} pairs, each q the LITERAL ntVec
   momentum (an env key by construction: momentumOf collected every ntVec momentum into the env basis),
   so legStr maps it straight to a vlc without re-decomposing over a non-atomic basis. A μ carried by
   no ntVec ⇒ a free (open) gluon id {"free",μ}; an ill-formed leg (a μ-bearing factor that is not a
   plain coeff·ntVec sum) ⇒ $Failed (stay conservative — the commutator then distributes as before). *)
ntSigmaLeg[leg_, termFactors_] := Module[{mu = First[leg], momFactors, terms},
  momFactors = Times @@ Select[termFactors, ! FreeQ[#, ntVec[_, mu]] &];
  If[momFactors === 1, Return[{"free", mu}]];
  terms = If[Head[momFactors] === Plus, List @@ momFactors, {momFactors}];
  Catch[
    {"slash", Sort[Function[trm, Module[{q, c},
       q = Cases[trm, ntVec[qq_, m_] /; m === mu :> qq, {0, Infinity}];
       c = trm /. ntVec[_, m_] /; m === mu :> 1;
       If[Length[q] =!= 1 || ! FreeQ[c, ntVec] || ! NumericQ[c], Throw[$Failed, "sigleg"]];
       {c, First[q]}]] /@ terms]},
    "sigleg"]];

(* chain-order the 2 gammas of one commutator term (spinor in->out) → {legFirst, legSecond, din, dout,
   scalar}, or $Failed if the term is not a clean isolated 2-γ pair. *)
ntSigmaTermInfo[term_] := Module[{tf, gs, g1, g2, first, second, din, dout},
  tf = If[Head[term] === Times, List @@ term, {term}];
  gs = Cases[tf, _ntGamma];
  If[Length[gs] =!= 2, Return[$Failed]];
  (* the two γ's (and their ntVecs) must be the term's ONLY tensor structure — no colour, projector,
     metric, other Dirac heads: anything else means the Plus is not a clean bare commutator. *)
  If[! FreeQ[tf, _ntGamma5 | _ntDeltaDirac | _ntSigma | _ntSUNT | _ntSUNDeltaFund | _ntSUNf | _ntSUNDeltaAdj |
                _ntTransProj | _ntLongProj | _ntMetric | _ntEpsilon],
    Return[$Failed]];
  {g1, g2} = gs;
  Which[
    g1[[3]] === g2[[2]], first = g1; second = g2; din = g1[[2]]; dout = g2[[3]],
    g2[[3]] === g1[[2]], first = g2; second = g1; din = g2[[2]]; dout = g1[[3]],
    True, Return[$Failed]];
  {ntSigmaLeg[first, tf], ntSigmaLeg[second, tf], din, dout, Times @@ Select[tf, scalarQ]}];

(* {scalar, ntSigma[...]} for a 2-term Plus that is a bare γ-commutator, else $Failed. The two terms
   must share spinor endpoints, carry the SAME two legs in SWAPPED order, and have opposite scalar.
   Result uses term-1's leg order and scalar: c·([X,Y] block) (sign-symmetric — either term gives the
   same fold, since c_2 = −c_1 and swapping legs negates the bracket). *)
ntRecognizeComm[p_] := Module[{terms, a, b},
  If[Head[p] =!= Plus || Length[p] =!= 2, Return[$Failed]];
  terms = List @@ p;
  a = ntSigmaTermInfo[terms[[1]]]; b = ntSigmaTermInfo[terms[[2]]];
  If[a === $Failed || b === $Failed, Return[$Failed]];
  If[a[[3]] =!= b[[3]] || a[[4]] =!= b[[4]], Return[$Failed]];      (* same spinor endpoints *)
  If[! (a[[1]] === b[[2]] && a[[2]] === b[[1]]), Return[$Failed]];  (* legs swapped between terms *)
  If[a[[1]] === a[[2]], Return[$Failed]];                          (* identical legs ⇒ [X,X]=0 *)
  If[Simplify[a[[5]] + b[[5]]] =!= 0, Return[$Failed]];            (* opposite scalar sign *)
  {a[[5]], ntSigma[a[[1]], a[[2]], a[[3]], a[[4]]]}];

(* replace every recognized bare γ-commutator Plus factor with its single ntSigma token + scalar.
   Set env NT_NO_SIGMA_FOLD=1 to disable the optimization (the commutator then distributes into two
   traces as before) — a safety switch and the baseline for measuring the fold's sub-term reduction. *)
$ntSigmaFold = (Environment["NT_NO_SIGMA_FOLD"] === $Failed);
(* Fold every bare γ-commutator Plus among `factors` into a single ntSigma token (gated by
   $ntSigmaFold), recursively: find a factor that is a Plus recognisable as a commutator [A,B]
   (ntRecognizeComm), replace it with its two σ legs, recurse until none remain — so the antisymmetric
   γ-pair is never distributed into two separate Dirac traces. *)
foldDiracSigma[factors_List] := Module[{commPlus, recognized},
  If[! $ntSigmaFold, Return[factors]];
  commPlus = FirstCase[factors, q_ /; (Head[q] === Plus && ntRecognizeComm[q] =!= $Failed), Missing[], {1}];
  If[MissingQ[commPlus], Return[factors]];
  recognized = ntRecognizeComm[commPlus];
  foldDiracSigma[Join[DeleteCases[factors, commPlus, {1}, 1], {recognized[[1]], recognized[[2]]}]]];

(* Split a colour-ENTANGLED Lorentz/Dirac component (the AqbqDirect147 4/7 quark-gluon vertex, where
   colour and Dirac do NOT factor globally, only per term) into per-colour-structure groups so the
   Lorentz-only builderInv can lower each. Algorithm:
     1. fold σ commutators, then EXPAND only the entangled Pluses (those mixing colour with
        Dirac/Lorentz) into branches — single-sector sums stay eager (no monomial blow-up).
     2. each branch → {colourProduct, {core, scal, restStr}} via compileDirac / chunkLorInv.
     3. GatherBy colour product, emitting ONE {colourNet, bodyNets, scalar, restNets} entry per group
        (so the net count tracks the colour graph, not the branch×ordering explosion); the generator
        sums each group's branches at runtime. *)
splitColourGroupsInv[factors0_, ids_, env_, mask_, nc_] := Module[
  {factors = foldDiracSigma[factors0],
   entangledQ, needExpand, keepAll, distributed, terms, branchNets, groups},
  entangledQ[x_] := Head[x] === Plus && (colourEntangledQ[x] || ! FreeQ[x, _ntGamma | _ntGamma5 | _ntDeltaDirac]);
  needExpand = Select[factors, entangledQ];
  keepAll    = Select[factors, ! entangledQ[#] &];
  distributed = Expand[Times @@ needExpand];          (* small: product of the entangled Pluses only *)
  terms = If[Head[distributed] === Plus, List @@ distributed, {distributed}];
  (* per branch -> {colourProduct, list-of-{bodyOrCore,scal,restStr}} *)
  branchNets = Function[term, Module[{termFactors, colProd, rest},
     termFactors = Join[If[Head[term] === Times, List @@ term, {term}], keepAll];
     colProd = Times @@ Cases[termFactors, Alternatives @@ ctHeadsInv];
     rest = DeleteCases[termFactors, Alternatives @@ ctHeadsInv];   (* gammas + Lorentz + numeric coeff (no colour) *)
     {colProd, If[! FreeQ[rest, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntDressedNum],
        {compileDirac[rest, ids, env, mask, nc]},          (* gamma chain: {core, scal, projectorRest} *)
        ({#[[1]], #[[2]], ""} &) /@ chunkLorInv[Times @@ rest, ids, env, mask, nc]]}]] /@ terms;
  groups = GatherBy[branchNets, First];
  (* one {colourNet, bodyNet, scalar, restNet} entry per colour group.
     INV backend: fuse a group's branches into ONE NetVal body via add/sc/contract — valid because
     each branch's `core` is itself a NetVal (dirac_value(...)); restNet is a string.
     NUMERIC backend: a branch's `core` is a BARE `DiracNet{...}` (contracted against its Lorentz
     rest only at runtime by numeric_value_netval), so fusing it into `add(sc(DiracNet{...}), …)` is
     a C++ type error (the `1/4/7` bug). Keep ONE entry per colour group (so the net count tracks the
     graph, not the branch×ordering explosion), but carry the group's branches as a LIST of sub-terms:
     bodyNet = {core_b…} (each a DiracNet literal, or a Lorentz NetVal for a gamma-free branch),
     restNet = {{rest_b, scal_b}…} (parallel). The generator sums mp[i]=Σ_b scal_b·numeric_value_netval(
     dnet_b, lnet_b) — matching how the inv body sums its branches, but type-correct for numeric. *)
  Function[group, Module[{colProd = group[[1, 1]], colNet, colScalar, branchRecs},
     {colNet, colScalar} = If[colProd === 1, {"SUNNet{}", 1}, compileColG[colProd, ids]];
     branchRecs = Flatten[group[[All, 2]], 1];                   (* each = {core, scal, restStr} *)
     {colNet, branchRecs[[All, 1]], colScalar, ({#[[3]], #[[2]]} &) /@ branchRecs}]] /@ groups];

(* ---- colour/group dialect: a constant SU(N) component is a product of structure constants,
        generators, and Kronecker deltas; emit it as a `SUNNet` literal for the generator's
        numeric SU(N) contraction (the et engine can't instantiate the four-gluon colour tensor
        type). Each head carries its own group rank N as the leading argument, so one net can mix
        several groups (e.g. colour SU(Nc) ⊗ flavour SU(Nf)) — sun_net.hpp's sun_value_cx
        contracts each rank separately and multiplies. Returns {colourNetString, factoredScalar}. *)
colFacG[ntSUNf[n_, a_, b_, c_], ids_] :=
  "SUN::f(" <> ToString[n] <> "," <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> "," <> ToString[ids[c]] <> ")";
colFacG[ntSUNDeltaAdj[n_, a_, b_], ids_] :=
  "SUN::deltaAdj(" <> ToString[n] <> "," <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> ")";
colFacG[ntSUNT[n_, a_, i_, j_], ids_]      := "SUN::T(" <> ToString[n] <> "," <> ToString[ids[a]] <> "," <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> ")";
colFacG[ntSUNDeltaFund[n_, i_, j_], ids_]  := "SUN::deltaFund(" <> ToString[n] <> "," <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> ")";
(* per-component diagonal dressings: register the (name, scale) -> dressing id, emit a diag factor. *)
colFacG[ntSUNDiagFund[n_, i_, j_, name_, scale_], ids_] :=
  "SUN::diagFund(" <> ToString[n] <> "," <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> "," <> ToString[diagDrId[name, scale, n, "fund"]] <> ")";
colFacG[ntSUNDiagAdj[n_, a_, b_, name_, scale_], ids_] :=
  "SUN::diagAdj(" <> ToString[n] <> "," <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> "," <> ToString[diagDrId[name, scale, n^2 - 1, "adj"]] <> ")";
compileColG[e_, ids_] := Module[{parts = If[Head[e] === Times, List @@ e, {e}], sc, tn},
  (* a colour/flavour factor raised to an integer power (e.g. deltaAdjFlav^2 from a CLOSED meson
     loop: delta_adj(a,b)^2 -> the flavour trace N^2-1) must be expanded into repeated SUNNet
     factors so the C++ sun_value contracts the shared indices — colFacG handles a single head,
     not Power[head,k], so an un-expanded power leaks Mathematica syntax into the generator. *)
  parts = parts /. Power[b_, k_Integer?Positive] /; ! scalarQ[b] :> Sequence @@ ConstantArray[b, k];
  sc = Select[parts, scalarQ]; tn = Select[parts, ! scalarQ[#] &];
  {"SUNNet{" <> StringRiffle[colFacG[#, ids] & /@ tn, ", "] <> "}", Times @@ sc}];

(* ---- et-vs-numeric size guard for a FUNDAMENTAL colour/flavour component ----------
   A small generator contraction (a quark loop's T^a) instantiates fine as an et ETensor and folds
   to a constexpr; a large one OOMs the compiler exactly like the four-gluon adjoint type. The metric
   below is the widest intermediate of contract_all under the greedy order — the product of the
   still-open axis extents (adjoint N^2-1, fundamental N) at each step — which directly predicts the
   ETensor type size. Above the threshold the component takes the numeric SUNNet path instead. *)
$NumTracerFundETMaxEntries = 200000;
labelDimAssoc[facs_] := Module[{assoc = <||>},
  Function[h, With[{n = sunRankOf[h]}, Switch[h,
      _ntSUNf | _ntSUNDeltaAdj | _ntSUNDiagAdj, (assoc[#] = n^2 - 1) & /@ labelsOf[h],
      _ntSUNT, (assoc[First[labelsOf[h]]] = n^2 - 1; (assoc[#] = n) & /@ Rest[labelsOf[h]]),
      _ntSUNDeltaFund | _ntSUNDiagFund, (assoc[#] = n) & /@ labelsOf[h]]]] /@ Select[facs, tensorQ];
  assoc];
fundMetric[facs_] := Module[
  {tn = orderFactors[Select[facs, tensorQ]], dimA = labelDimAssoc[facs], cnt = <||>, widest = 1, open},
  Do[(cnt[#] = Lookup[cnt, #, 0] + 1) & /@ labelsOf[f];
     open = Keys[Select[cnt, OddQ]];
     widest = Max[widest, Times @@ (dimA /@ open)], {f, tn}];
  widest];

(* ---- Dirac trace -> Lorentz NetVal (for the invariant-basis path) ----------
   The et engine takes a spinor trace by closing a cyclic chain; the inv generator works in
   scalar products, so we instead apply the gamma-trace THEOREM symbolically, turning the
   closed chain of gammas into a signed sum of metric products. Each slashed momentum is
   ntGamma[mu,..] * ntVec[q,mu], so the metric mu-pairings contract the vecs/projectors that
   share those Lorentz labels — leaving a pure-Lorentz network the existing reduce->rebase
   pipeline handles unchanged. (No gamma5 support yet; gamma5-bearing traces are not expanded.) *)
diracIn[ntGamma[_, a_, _]]      := a;  diracOut[ntGamma[_, _, b_]]      := b;
diracIn[ntGamma5[a_, _]]        := a;  diracOut[ntGamma5[_, b_]]        := b;
diracIn[ntSigma[_, _, a_, _]]   := a;  diracOut[ntSigma[_, _, _, b_]]   := b;
diracIn[ntDeltaDirac[a_, _]]    := a;  diracOut[ntDeltaDirac[_, b_]]    := b;
diracIn[ntDressedNum[_, a_, _]] := a;  diracOut[ntDressedNum[_, _, b_]] := b;

(* walk the closed spinor loop, returning trace-ordered tokens: a γ's Lorentz label, or the marker
   $g5 for a γ5 (a spinor-δ just connects indices and contributes no token). *)
orderDiracChain[facs_] := Module[{byIn, cur, toks = {}, seen = {}},
  byIn = Association[(diracIn[#] -> #) & /@ facs];
  cur = First[facs];
  While[! MemberQ[seen, cur],
    AppendTo[seen, cur];
    Which[MatchQ[cur, _ntGamma], AppendTo[toks, First[cur]],
          MatchQ[cur, _ntGamma5], AppendTo[toks, $g5]];
    cur = Lookup[byIn, diracOut[cur], Missing[]];
    If[MissingQ[cur], Break[]]];
  toks];

(* tr[g^m1 ... g^m2n] / 4 = sum over pairings of signed products of metrics (recursive). The pairing
   STRUCTURE depends only on the chain LENGTH, not the labels: gammaTraceSum[{a,b,c,d}] is
   gammaTraceSum[{1,2,3,4}] with positions -> labels. A quark triangle distributes into hundreds of
   branches that each trace a long (10-12 gamma) chain — recomputing the (2n-1)!! pairing sum per
   branch was catastrophic (minutes / GBs per component). Instead compute the canonical positional
   trace ONCE per length (memoized) and substitute each branch's labels — hundreds of cheap
   substitutions replace hundreds of exponential recomputations. *)
rawGammaTrace[{}] := 1;
rawGammaTrace[lst_] := Sum[(-1)^j ntMetric[lst[[1]], lst[[j]]] rawGammaTrace[Delete[lst, {{1}, {j}}]],
  {j, 2, Length[lst]}];
canonGammaTrace[n_] := canonGammaTrace[n] = rawGammaTrace[Range[n]];   (* positions 1..n, memoized per length *)
gammaTraceSum[{}] := 1;
gammaTraceSum[lst_] := canonGammaTrace[Length[lst]] /. Dispatch[Thread[Range[Length[lst]] -> lst]];

(* tr[γ5 g^m1 ... g^m2n] / 4 = the chiral trace. Vanishes for < 4 gammas; the 4-gamma base is the
   Levi-Civita ε — the normalisation tr(γ5 γ γ γ γ) = 4 ε is pinned against the typed-out component
   engine (test_inv_eps, constant C = 4). ≥ 6 gammas peel one metric exactly like the ordinary trace,
   with the ε base seeding the recursion (the standard FORM chiral reduction). *)
gammaTraceSum5[lst_] := Which[
  Length[lst] < 4, 0,
  Length[lst] == 4, ntEpsilon[lst[[1]], lst[[2]], lst[[3]], lst[[4]]],
  True, Sum[(-1)^j ntMetric[lst[[1]], lst[[j]]] gammaTraceSum5[Delete[lst, {{1}, {j}}]], {j, 2, Length[lst]}]];

(* trace of a token chain: collapse any number of γ5 to ≤ 1 (γ5 γ^μ = −γ^μ γ5, γ5² = 1). The sign to
   slide every γ5 together past the intervening γ's is (−1)^(#{γ before each γ5}); an even γ5 count →
   the ordinary trace, an odd count → the chiral trace. (This is what guarantees the C++ ε-network's
   "≤ 1 ε per term" invariant.) *)
diracTrace[toks_] := Module[{n5 = Count[toks, $g5], gammas = DeleteCases[toks, $g5], sign},
  sign = (-1)^Total[(Length[Cases[Take[toks, # - 1], Except[$g5]]]) & /@ Flatten[Position[toks, $g5]]];
  sign * If[EvenQ[n5], gammaTraceSum[gammas], gammaTraceSum5[gammas]]];

(* expand a Dirac-trace component (factor list) to a pure-Lorentz expression; non-Dirac
   components pass through unchanged. The tr(1)=4 is the explicit leading factor.
   LEGACY: the inv backend now traces in C++ (compileDirac → network::dirac_value); this symbolic
   path is kept for the ET/Dense backends and for cross-validation. *)
expandDiracComponent[factors_List] := Module[{dir, rest},
  dir = Select[factors, MatchQ[#, _ntGamma | _ntGamma5 | _ntDeltaDirac] &];
  rest = DeleteCases[factors, _ntGamma | _ntGamma5 | _ntDeltaDirac];
  If[dir === {}, Times @@ factors,
    4 diracTrace[orderDiracChain[dir]] (Times @@ rest)]];

(* ---- Dirac trace in the C++ generator (et/inv/dirac.hpp `dirac_value`) -----------------------
   The inv backend emits the gamma chain as a `DiracNet` literal and the generator contracts the
   closed spinor loop NUMERICALLY into a Lorentz net (metrics/vectors over the free gluon legs), then
   `reduce` folds it with the projector — instead of expanding the (2n−1)!! pairing sum symbolically
   in Mathematica (which exploded for the 3-gluon-vertex quark triangle). This parallels how colour is
   folded by `SUNNet`/`compileColG`. *)
$ntDiracFree = 900000; (* fresh Lorentz-label base for slash–slash pairings: above every component id *)

(* the gamma/gamma5 factors in trace order (walk the closed spinor loop; spinor-δ connectors carry no
   token but are followed). Like `orderDiracChain` but returns the factors so each can be classified.

   UNDIRECTED cycle walk: each Dirac factor is an EDGE between its two spinor labels (spinorLabelsHead);
   every label in a closed spinor loop has degree 2, so the walk is deterministic. We seed at
   `First[facs]` ENTERING on its `diracIn` label (so it leaves on `diracOut`) and at each step leave the
   current factor by its OTHER endpoint, picking the unique unvisited neighbour there. For a
   consistently-oriented loop (every node one in / one out — Zq/ZA/ZA3/ZA4/ZAqbq, which close via an
   oriented γ) this reproduces the old directed `byIn` forward walk token-for-token (byte-identical).
   The robustness is for a loop closed by a SYMMETRIC external spinor-δ (`ntDeltaDirac[d1,d2]`, e.g. the
   σL scalar external projector): the old directed walk collided on the shared `diracIn`/`diracOut` keys
   and could not traverse the cycle, collapsing the trace to a single γ5 → tr(γ5)=0. The undirected
   adjacency closes the cycle correctly regardless of start/orientation (the trace is cyclic). *)
orderDiracFacs[facs_] := Module[{nodeFacs = Association[], cur = 1, prevLabel, out = {}, seen = {},
    labels, exitLabel, nexts, nTok},
  Do[Module[{ls = spinorLabelsHead[facs[[i]]]},
       (nodeFacs[#] = Append[Lookup[nodeFacs, #, {}], i]) & /@ ls], {i, Length[facs]}];
  prevLabel = diracIn[facs[[1]]];   (* enter First on its diracIn so it exits on diracOut *)
  While[! MemberQ[seen, cur],
    AppendTo[seen, cur];
    If[MatchQ[facs[[cur]], _ntGamma | _ntGamma5 | _ntSigma | _ntDressedNum], AppendTo[out, facs[[cur]]]];
    labels = spinorLabelsHead[facs[[cur]]];
    exitLabel = First[DeleteCases[labels, prevLabel], Missing[]];   (* the OTHER endpoint *)
    If[MissingQ[exitLabel], Break[]];
    nexts = Select[DeleteCases[Lookup[nodeFacs, exitLabel, {}], cur], ! MemberQ[seen, #] &];
    If[nexts === {}, Break[]];
    prevLabel = exitLabel; cur = First[nexts]];
  (* LOUD GUARD: a closed spinor loop must consume EVERY token-bearing factor. If the walk fell short
     (a fragmented / improperly-closed loop), the emitted trace would silently drop γ structure (the
     historical hSigL meson-sector collapse: a single surviving γ5 → tr(γ5)=0). Abort rather than emit a
     wrong kernel. The δ-only "connector" factors carry no token, so they are excluded from the count. *)
  nTok = Count[facs, _ntGamma | _ntGamma5 | _ntSigma | _ntDressedNum];
  If[Length[out] =!= nTok,
    Print["[NumTracer] ERROR: orderDiracFacs walk consumed ", Length[out], " of ", nTok,
          " token-bearing Dirac factors — a spinor loop did not close (would silently drop γ ",
          "structure → collapsed trace). Loop factors:\n  ", facs];
    Abort[]];
  out];

(* A component may contain SEVERAL independent closed spinor loops (a quark loop + the
   projection-closed external line, tied together only by gluon propagators). orderDiracFacs walks ONE
   loop and would drop the rest; partition the Dirac factors into spinor-connected groups first
   (two factors share a loop iff they share a spinor index), then order each loop. Returns a LIST of
   ordered token lists (one per loop) — a single-loop component yields a one-element list. The C++
   contraction traces each loop separately and contracts their shared gluon legs (DFac::LoopSep). *)
orderDiracLoops[facs_] := If[Length[facs] <= 1, {orderDiracFacs[facs]},
  Module[{sp = spinorLabelsHead /@ facs, edges, g, comps},
    edges = Select[Subsets[Range[Length[facs]], {2}], IntersectingQ[sp[[#[[1]]]], sp[[#[[2]]]]] &];
    g = Graph[Range[Length[facs]], UndirectedEdge @@@ edges];
    comps = ConnectedComponents[g];
    (* a SINGLE spinor loop → the exact original path (orderDiracFacs[facs]), so existing flows stay
       byte-identical; only genuinely multi-loop components are split (indices sorted to preserve the
       original relative order within each loop). *)
    If[Length[comps] <= 1, {orderDiracFacs[facs]},
      (orderDiracFacs[facs[[Sort[#]]]]) & /@ comps]]];

(* {netString, scalar} for a colour-free component (gammas + Lorentz + numeric coeff). The gamma chain
   becomes `dirac_value(DiracNet{...})` (each γ: a SLASH `dslash` if its Lorentz index is contracted
   with an `ntVec[q,μ]`, else a FREE leg `dgamma(μ)` that contracts the projector); the remaining
   (Lorentz) factors compile through `compileTInv` and are `contract`ed with the trace. A component
   with no Dirac factor is just its Lorentz net. *)
(* frame resolver for dressed-numerator option coefficients (ntSP/ntVec[q,i] -> components). Set in
   mkGenerateKernel to the diagram's resolveScale; Identity when the dressed path is inactive. *)
$ntDressResolve = Identity;
(* the C++ DSlot literal for one ntDressedNum: a sum of DSlotOpt{Cx coeff, {dress ids}, slash?, vlc}.
   Each option's scalar coefficient is frame-resolved then split into a complex number × dressing atoms
   (drDecompose); a "slash" structure's momenta become a vlc of {1.0, env Base} pairs (coeff-1, like the
   single-vec slash above), an "ident" structure is the spinor identity (slash=false, empty vlc). *)
dressedSlotStr[ntDressedNum[opts_, _, _], env_] := "DSlot{" <> StringRiffle[
  Function[opt, Module[{num, dr, vlcStr},
     {num, dr} = drDecompose[$ntDressResolve[opt[[1]]]];
     vlcStr = If[opt[[2, 1]] === "slash",
        StringRiffle[Function[cq, Module[{c = cq[[1]], q = cq[[2]]},
           If[! KeyExistsQ[env, q], Print["[NumTracer] ERROR: dressed slash momentum ", q,
              " absent from env ", Keys[env]]; Abort[]];
           "{" <> cppNum[c] <> "," <> ToString[env[q]["Base"]] <> "}"]] /@ opt[[2, 2]], ", "], ""];
     "DSlotOpt{Cx{" <> cppNum[Re[num]] <> "," <> cppNum[Im[num]] <> "}, {" <>
       StringRiffle[ToString /@ dr, ", "] <> "}, " <>
       If[opt[[2, 1]] === "slash", "true", "false"] <> ", {" <> vlcStr <> "}}"]] /@ opts, ", "] <> "}";

(* Turn a component's factor list into the C++ Dirac-chain token(s) + its Lorentz "rest".
   Algorithm:
     1. pull out the Dirac heads (γ/γ5/σ/δ/dressed-numerator); if none, defer to compileTInv.
     2. order them into independent spinor loops (orderDiracLoops) and map each to a token string
        (tokenOf): a γ leg becomes dgamma/dslash, a σ becomes the matching dcomm* builder (legStr),
        a dressed-numerator becomes a SLOT (dtslot) whose DSlot is recorded.
     3. LOUD-GUARD that no Dirac head leaked into the non-Dirac "rest" (an un-collected dressed sum) —
        that would silently collapse/leak the trace; abort if so.
     4. compile the Lorentz rest separately (restCompiled = {restStr, scalar}) and return the Dirac
        core, scalar and rest SEPARATELY so the contraction stays deferred (rest emitted once, shared
        across a colour group). A dressed chain returns the ntDressedCore[chain, slots] marker instead. *)
compileDirac[factors_, ids_, env_, mask_, nc_] := Module[
  {diracFacs, loops, loopStrs, slashVecs = {}, tokenOf, restFacs, restCompiled, legStr, sigStr, slots = {}, slotN = 0, dressed},
  diracFacs = Select[factors, MatchQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntDressedNum] &];
  If[diracFacs === {}, Return[Append[compileTInv[Times @@ factors, ids, env, mask, nc], ""]]];
  loops = orderDiracLoops[diracFacs];   (* one ordered token list per independent spinor loop *)
  dressed = ! FreeQ[diracFacs, _ntDressedNum];
  (* an ntSigma leg → its C++ arg string. A slashed leg is a list of {coeff, q} pairs (q a literal
     ntVec momentum, hence an env key) ⇒ a vlc `{{c1,b1},{c2,b2},…}` of (coeff, env Base) pairs (a
     single atomic momentum reduces to `{{1.,Base}}`); a free leg is its open Lorentz id. Guard loudly
     if a leg momentum is somehow not in env (it always should be — momentumOf collected it). *)
  legStr[{"slash", pairs_List}] := "{" <> StringRiffle[
     Function[cm, Module[{c = cm[[1]], q = cm[[2]]},
        If[! KeyExistsQ[env, q],
          Print["[NumTracer] ERROR: σ slash leg momentum ", q, " absent from env ", Keys[env]]; Abort[]];
        "{" <> cppNum[c] <> "," <> ToString[env[q]["Base"]] <> "}"]] /@ pairs, ", "] <> "}";
  legStr[{"free", mu_}]   := ToString[ids[mu]];
  (* the bare-commutator [A,B] builder, picked by each leg's free/slash kind (see et/inv/dirac.hpp). *)
  sigStr[legA_, legB_] := Module[{ta = legA[[1]], tb = legB[[1]], sa = legStr[legA], sb = legStr[legB]},
    Which[
      ta === "free"  && tb === "free",  "dcomm("    <> sa <> ", " <> sb <> ")",
      ta === "slash" && tb === "slash", "dcomm_ss(" <> sa <> ", " <> sb <> ")",
      ta === "free"  && tb === "slash", "dcomm_fs(" <> sa <> ", " <> sb <> ")",
      True,                             "dcomm_sf(" <> sa <> ", " <> sb <> ")"]];
  (* one chain token. A ntDressedNum becomes a SLOT (dtslot index, its DSlot recorded); the others are
     the usual gamma/slash/γ5/σ tokens. In the dressed case every fixed token is wrapped in dtfix(...)
     so the chain is a std::vector<DChainTok>; otherwise the bare token goes into a DiracNet{...}. *)
  tokenOf = Function[gf,
     Which[
       MatchQ[gf, _ntDressedNum],
         (AppendTo[slots, dressedSlotStr[gf, env]];
          With[{k = slotN}, slotN++; "dtslot(" <> ToString[k] <> ")"]),
       MatchQ[gf, _ntGamma5], If[dressed, "dtfix(dg5())", "dg5()"],
       MatchQ[gf, _ntSigma],  With[{s = sigStr[gf[[1]], gf[[2]]]}, If[dressed, "dtfix(" <> s <> ")", s]],
       True,
       Module[{mu = First[gf], vq, t},
         vq = Cases[factors, ntVec[q_, m_] /; m === mu :> q];
         t = If[vq =!= {},
           (AppendTo[slashVecs, ntVec[First[vq], mu]];
            "dslash({{1.0," <> ToString[env[First[vq]]["Base"]] <> "}})"),
           "dgamma(" <> ToString[ids[mu]] <> ")"];
         If[dressed, "dtfix(" <> t <> ")", t]]]];
  (* one token-string per spinor loop (mapped in order so the slot/slashVec side effects accumulate) *)
  loopStrs = (StringRiffle[tokenOf /@ #, ", "]) & /@ loops;
  restFacs = DeleteCases[factors, Alternatives @@ Join[diracFacs, slashVecs]];
  (* LOUD GUARD against the silent dressed-numerator fall-through. `diracFacs` captures only BARE Dirac
     heads; a dressed propagator-numerator sum (Mq·δ − I·Z·γ·p̸) that was NEITHER distributed
     (expandBridges) NOR collected (rewriteDressedNums → ntDressedNum) survives as a `Plus`/`Power`
     factor with Dirac heads NESTED inside it. Such a factor lands in `restFacs` → `compileTInv`, which
     handles only Lorentz/colour/scalars: it then SILENTLY drops the nested γ/δ (the trace collapses,
     e.g. a pion loop → tr(γ5)=0) or leaks the raw head into the generated C++ (undeclared `ntGamma` /
     `ntDeltaDirac`). Either way the diagram is wrong with no warning (observed in the hSigL meson
     sector, ~1% of FORM). Detect it and abort with the offending factor rather than emit a wrong kernel.
     A correctly-handled diagram has all Dirac structure in `diracFacs` (bare heads or ntDressedNum), so
     `restFacs` is Dirac-free; this never trips on the validated flows. *)
  With[{leak = Select[restFacs, ! FreeQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntDressedNum] &]},
    If[leak =!= {},
      Print["[NumTracer] ERROR: un-handled Dirac structure in a non-Dirac factor — a dressed ",
            "propagator-numerator sum was NEITHER distributed NOR collected into ntDressedNum, so the ",
            "numeric backend would silently drop/leak its gamma structure (collapsed trace or ",
            "untranslated C++). This is a front-end collection gap (collectibleDiracSumQ rejected a sum ",
            "that distributeQ also skipped). Offending factor(s):\n  ", leak];
      Abort[]]];
  restCompiled = compileTInv[Times @@ restFacs, ids, env, mask, nc];   (* {restStr, scal}; restStr "" if trivial *)
  (* return the Dirac core and the Lorentz "rest" (projector) SEPARATELY: {coreStr, scal, restStr}.
     The contract is deferred — splitColourGroupsInv factors a rest shared across a colour group out of
     the Dirac-trace sum so it is emitted ONCE (a shared sub-net) and contracted lazily in the
     generator, never materialising the |trace|×|projector| product. restStr=="" ⇒
     no projector. The numeric backend emits the bare DiracNet{...}; the runtime contracts it against
     the Lorentz rest via numeric_value_netval (the matrix-product trace). A DRESSED chain (one or more
     ntDressedNum slots) is returned as the marker ntDressedCore[chainStr, slotsStr]; the generator
     contracts it via numeric_value_dressed_netval (one collected DPoly trace). *)
  (* join the per-loop token strings with a loop separator so the runtime traces each spinor loop
     independently and contracts their shared gluon legs; a single-loop component has no separator and
     emits the identical net as before. *)
  If[dressed,
    Module[{chain = "std::vector<DChainTok>{" <> StringRiffle[loopStrs, ", dtfix(dloopsep()), "] <> "}",
            slotV = "std::vector<DSlot>{" <> StringRiffle[slots, ", "] <> "}"},
      {ntDressedCore[chain, slotV], restCompiled[[2]], restCompiled[[1]]}],
    Module[{core = "DiracNet{" <> StringRiffle[loopStrs, ", dloopsep(), "] <> "}"},
      {core, restCompiled[[2]], restCompiled[[1]]}]]];

(* ---- numeric (matrix-product) backend: component table + user symbols (task #22) -------------
   The numeric backend has NO sp-invariant basis: it evaluates scalar products numerically from each
   momentum's 4 COMPONENTS. So instead of an sp-invariant decomposition it needs
   only: the polynomial variables (the free user symbols), each fundamental momentum's 4 components as
   MPoly-builder C++, and the C++ fill formula for each symbol (a kernel arg, or a derived symbol like
   sin1 = sqrt(1-cos1^2)). Composite momenta resolve by component arithmetic via resolveComponents,
   exactly like the frame path — the user "Components" association is just a frame whose entries are
   polynomial. *)

(* Polynomialise a frame whose loop direction carries Sqrt[1-cos^2]: introduce a sin symbol per angle
   and record its definition. A user-supplied (already polynomial) "Components" assoc is a no-op. *)
(* Unit-loop spec (the COMPACT numeric parametrisation that matches the sp invariant count). The
   naive frame fallback bloats because the LOOP momentum is decomposed into angle PRODUCTS
   (l1·cos1, l1·sin1·cos2, … — degree-3), so every l·p_i expands into many angle monomials and the
   degrees compound through a projector chain. Instead write the loop as MAGNITUDE × UNIT-DIRECTION:
   comp_μ = l1 · Uμ where Uμ are degree-1 symbols (fill computes Uμ = dirμ from the kernel angles).
   Then l·p_i = l1·p·Σ Uμ u_i^μ is degree-1 in Uμ (as compact as sp(l,p_i)), and the bare-loop
   denominator Σ(l1·Uμ)² = l1²·ΣUμ² collapses to the MONOMIAL l1² under the unit constraint ΣUμ²=1
   (returned as a unit `group`), so its 1/l² atom cancels — exactly like inv's `rel`. Externals depend
   only on `p` (kept as the numeric p-vector). The loop is the momentum whose components carry `magSym`. *)
unitLoopFrameSpec[frame_, pSym_, magSym_] := Module[{defs = <||>, groups = {}, n = 0, nf},
  nf = Association @ KeyValueMap[Function[{q, comps},
     q -> If[SubsetQ[{pSym}, Variables[comps]], comps,         (* external: numeric p-vector *)
        Module[{grp = {}, nc},                                  (* loop: comp_μ = magSym · Uμ *)
          nc = Table[Module[{dir = Coefficient[comps[[mu]], magSym], s},
                If[dir === 0, 0,
                   s = Symbol["ntU$" <> ToString[n++]]; defs[s] = dir; AppendTo[grp, s]; magSym s]], {mu, 1, 4}];
          AppendTo[groups, grp]; nc]]], frame];
  {nf, defs, groups}];

(* General frame -> {polyFrame, defs}: replace the trig sub-expressions in the frame components with
   fresh symbols so the components become polynomials in (magnitudes, cos's, sin's). Introduces a
   symbol for each loop polar factor Sqrt[1-cos^2] (a sin) and each bare-angle Cos/Sin (the 4-point
   frame's φ); `defs` maps each new symbol back to its closed form (so the kernel can fill it). The
   fallback path when the unit-loop parametrisation (unitLoopFrameSpec) does not apply. *)
polyFrameSpec[frame_] := Module[{defs = <||>, vals, rules = {}, polyFrame},
  vals = PowerExpand[Values[frame]];
  (* loop polar factor Sqrt[1-cos^2] -> a sin symbol *)
  Do[With[{s = Symbol["ntSin$" <> SymbolName[c]]}, defs[s] = Sqrt[1 - c^2]; AppendTo[rules, Sqrt[1 - c^2] -> s]],
    {c, DeleteDuplicates @ Cases[vals, Power[1 - c_^2, 1/2] :> c, Infinity]}];
  (* trig of a bare angle symbol (e.g. Cos[phi], Sin[phi] in the 4-point frame) -> cos/sin symbols *)
  Do[With[{s = Symbol["ntCos$" <> SymbolName[a]]}, defs[s] = Cos[a]; AppendTo[rules, Cos[a] -> s]],
    {a, DeleteDuplicates @ Cases[vals, Cos[a_Symbol] :> a, Infinity]}];
  Do[With[{s = Symbol["ntSinA$" <> SymbolName[a]]}, defs[s] = Sin[a]; AppendTo[rules, Sin[a] -> s]],
    {a, DeleteDuplicates @ Cases[vals, Sin[a_Symbol] :> a, Infinity]}];
  polyFrame = Association @ Thread[Keys[frame] -> (vals /. rules)];
  {polyFrame, defs}];

(* Whether `frame` matches the unit-loop spec's assumption: every momentum is either an external
   depending only on `pSym`, or a loop whose every component is `dir·magSym` (proportional to the
   single magnitude). A finite-T frame breaks this — the external carries an independent temporal
   component p0 and the loop an independent l0 (neither ∝ magSym) — so we must fall back to the
   general polyFrameSpec. *)
unitLoopOkQ[frame_, pSym_, magSym_] := AllTrue[Values[frame], Function[comps,
  Module[{cc = PowerExpand[comps]},
    SubsetQ[{pSym}, Variables[cc]] ||
    AllTrue[Range[4], (Simplify[cc[[#]] - Coefficient[cc[[#]], magSym] magSym] === 0) &]]]];

numericComponents[env_, frame_, symDefs_, unitGroups_:{}] := Module[
  {compExpr, usyms, nsym, mpcpp, compCpp, varFill, vfill, idx, units},
  (* 4 components per momentum Base (polynomial in the user symbols). *)
  compExpr = Association @ KeyValueMap[#2["Base"] -> resolveComponents[#1, frame] &, env];
  (* reject only fractional powers of NON-numeric bases (e.g. Sqrt[1-cos1^2] not rewritten to a sin
     symbol); a numeric irrational coefficient like Sqrt[3]/2 from the 120-degree external frame is a
     perfectly good polynomial coefficient. *)
  With[{bad = DeleteDuplicates @ Cases[Values[compExpr], Power[b_, _Rational] /; ! NumericQ[b], Infinity]},
    If[bad =!= {},
      Print["numericComponents: non-polynomial components (fractional power of a symbol remains): ", bad];
      Abort[]]];
  usyms = Sort @ DeleteDuplicates @ Flatten[Variables /@ Values[compExpr]];
  nsym = Length[usyms];
  mpcpp[e_] := Module[{rules = CoefficientRules[e, usyms]},
    If[rules === {}, "MPoly(" <> ToString[nsym] <> ")",
      "(" <> StringRiffle[("MPoly::mono(" <> ToString[nsym] <> ",{" <>
          StringRiffle[ToString /@ #[[1]], ","] <> "},Cx{" <> cppNum[#[[2]]] <> ",0})") & /@ rules, " + "] <> ")"]];
  compCpp = Association @ KeyValueMap[#1 -> (mpcpp /@ #2) &, compExpr];
  vfill[s_] := If[KeyExistsQ[symDefs, s], cppFlat[symDefs[s]], SymbolName[s]];
  varFill = vfill /@ usyms;                          (* indexed by MPoly var id (0-based) *)
  (* Unit-constraint groups (ΣUμ²=1) as MPoly var-index lists, so the C++ reduce_units collapses the
     bare-loop denominator to the monomial l1² and the U·U projector factors. The caller supplies the
     groups as symbol lists (the loop's unit-direction components); drop any symbol not in usyms (a
     component that vanished) and any group with <2 surviving symbols. *)
  idx = AssociationThread[usyms -> Range[Length[usyms]] - 1];
  units = DeleteCases[Function[g, Lookup[idx, Select[g, KeyExistsQ[idx, #] &], Nothing]] /@ unitGroups,
                      _?(Length[#] < 2 &)];
  <|"nsym" -> nsym, "usyms" -> usyms, "compCpp" -> compCpp, "varFill" -> varFill,
    (* default -1 (-> comp size 0) for a purely scalar integrand whose 4-vector component env is
       empty: Max[{}] is -Infinity, which would leak into the C++ as the comp() vector size. The loop
       momentum of such a flow (e.g. a bosonic meson-potential tadpole) is carried by usyms/units. *)
    "symNamesCpp" -> (("(" <> # <> ")") & /@ varFill), "maxBase" -> Max[Append[#["Base"] & /@ Values[env], -1]],
    "units" -> units|>];

(* ---- emit the NUMERIC generator program. Builds each net's DiracNet chain + Lorentz/projector
        NetVal in C++, contracts them numerically (4×4 matrix products), folds colour per group, and
        PRINTS the committed straight-line kernel header. No reduce/rebase/ibp/sp-kinematics. *)
emitNumericGenerator[invNets_, invRest_, colourNets_, groups_, ncomp_, nsInner_, fillArgSig_,
    kns_:"numtracer_kernels", complexQ_:False, realOnlyG_:{}] := Module[
  {nNet = Length[invNets], nGrp = Length[groups], nsym = ncomp["nsym"], maxBase = ncomp["maxBase"],
   varFill = ncomp["varFill"], symNames = ncomp["symNamesCpp"], compCpp = ncomp["compCpp"], unitG = ncomp["units"],
   bb, tmpl, pre, unitPre, nUnits, units, decl, ddl, dStr, lStr, dscv, ddch, ddsl, hasDressed,
   dnetDefs, lnetDefs, dchDefs, dslDefs, allDefs, main, compInit, cseDefs, cseDecls},
  bb[x_] := ToString[x];
  (* DRESSED nets (symbolic dressing collection): a core may be ntDressedCore[chainStr, slotsStr]
     (a numerator structure-sum kept eager). hasDressed routes the generator to the DPoly branch:
     mp[] is a vector<DPoly>, a dressed sub-term contracts via numeric_value_dressed_netval, and the
     dressings ride the env as kind-2 `dress` leaves filled by fm.dress. The CSE / sub-term dedup
     (optimisations for the dense σ case) are SKIPPED when hasDressed so the non-dressed path stays
     byte-identical. *)
  hasDressed = ! FreeQ[invNets, _ntDressedCore];
  (* shared wrapper templates (same as the inv generator) so the net-builder strings compile. *)
  tmpl = "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal tproj(){ return projT(Mu,Nu,Lb,Inv); }\n" <>
    "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal lproj(){ return projL(Mu,Nu,Lb,Inv); }\n" <>
    "template<int Mu,int Nu,int Lb,int Mask,int InvS> NetVal mproj(){ return projM(Mu,Nu,Lb,InvS); }\n" <>
    "template<int Mu,int Nu,int Lb,int Mask,int Inv,int InvS> NetVal eproj(){ return projE(Mu,Nu,Lb,Inv,InvS); }\n" <>
    "template<int Mu,int Nu> NetVal lmetric(){ return met(Mu,Nu); }\n" <>
    "template<int Lbl,int Base,int Mask> NetVal lvec(){ return vec(Lbl,Base); }\n" <>
    "template<int A,int B,int C,int D> NetVal leps(){ return epsilon(A,B,C,D); }\n" <>
    "inline NetVal konst(double c){ return NetVal{PTerm{Cx{c,0}, {}}}; }\n" <>
    "template<class L> struct litco;\n" <>
    "template<numtracer::Cx C> struct litco<numtracer::expr::Lit<C>>{ static constexpr numtracer::Cx v=C; };\n" <>
    "template<class L> NetVal sc(NetVal x){ return scale(litco<L>::v, std::move(x)); }\n";
  (* per net: a colour group is a SUM of sub-terms. invNets[i] = {core_b…} (each a DiracNet literal
     for a gamma branch, or a Lorentz NetVal for a gamma-free branch); invRest[i] = {{rest_b,scal_b}…}
     parallel. Build per net the parallel lists of {DiracNet builder, NetVal builder} (a gamma-free
     branch → empty DiracNet + the whole net as the rest) plus the sub-term scalars; the generator
     sums mp[i] = Σ_b scal_b · numeric_value_netval(dn[i][b], ln[i][b]). *)
  {dStr, lStr, dscv, ddch, ddsl} = Transpose @ MapThread[Function[{cores, rss},
     If[cores === {}, {{}, {}, {}, {}, {}},
       Module[{ds, ls, dc, dl},
         {ds, ls, dc, dl} = Transpose @ MapThread[Function[{nv, rv},
            Which[
              MatchQ[nv, _ntDressedCore],                  (* dressed numerator: chain+slots, no DiracNet *)
                {"DiracNet{}", If[rv === "", "NetVal{}", rv], nv[[1]], nv[[2]]},
              StringMatchQ[nv, "DiracNet" ~~ ___],          (* gamma branch: DiracNet + projector rest *)
                {nv, If[rv === "", "NetVal{}", rv], "std::vector<DChainTok>{}", "std::vector<DSlot>{}"},
              True,                                         (* gamma-free branch: whole net is the rest *)
                {"DiracNet{}", nv, "std::vector<DChainTok>{}", "std::vector<DSlot>{}"}]],
            {cores, rss[[All, 1]]}];
         {ds, ls, rss[[All, 2]], dc, dl}]]], {invNets, invRest}];
  pre = StringJoin[
    "// GENERATED by MakeNTKernel — do not edit. Numeric matrix-product tensor traces.\n",
    "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n",
    "#include \"numtracer/codegen/gen.hpp\"\n#include \"numtracer/network/sun_net.hpp\"\n",
    "#include \"numtracer/numeric/numeric_contract.hpp\"\n#include \"numtracer/numeric/numeric_driver.hpp\"\n",
    "#include \"numtracer/expr/expr.hpp\"\n",
    "#include <iostream>\n#include <string>\n#include <utility>\n#include <vector>\n#include <array>\n#include <cstdlib>\n",
    "#include <thread>\n#include <atomic>\n#include <mutex>\n#include <chrono>\n#include <cstdio>\n#include <algorithm>\n#include <system_error>\n",
    "#include <sstream>\n#include <unordered_map>\n",
    "using numtracer::Cx;\n",
    "namespace numtracer::network {\n", tmpl, "}\n",
    "using namespace numtracer::network;\nusing namespace numtracer::numeric;\n"];
  unitPre = "// GENERATED by MakeNTKernel — do not edit. Numeric net-builder unit (compiled -O0).\n" <>
    "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n#include \"numtracer/expr/expr.hpp\"\n#include <utility>\n" <>
    (* dressed nets emit dch<i>()/dsl<i>() builders here (the big DChainTok/DSlot literals — moved OFF
       the single -O1 main TU onto these parallel -O0 units, since a 100k+-char braced-init is ~quadratic
       even at -O0); they need the dressed-token types from numeric_contract.hpp. Non-dressed units don't
       include it (stay byte-identical + fast). *)
    If[hasDressed, "#include \"numtracer/numeric/numeric_contract.hpp\"\n", ""] <>
    "using numtracer::Cx;\nnamespace numtracer::network {\n" <> tmpl <> "}\nusing namespace numtracer::network;\n" <>
    If[hasDressed, "using namespace numtracer::numeric;\n", ""];
  (* NET-LEVEL CSE: dense projections (e.g. the σ^μν struct-7 quark-gluon vertex) emit the SAME net
     sub-term thousands of times — measured 48624 lnet terms but only 66 distinct (534x), so the raw
     generator C++ balloons to ~25 MB and the -O0 compile dominates generation. Hash-cons each DISTINCT
     net term into a shared accessor `lc<k>()` / `dc<k>()` (a function-local `static const` so the net
     is also BUILT once at run time, not once per use), and reference it. Trivial/empty literals and
     unique terms stay inline. Correctness-preserving: each use copies the shared NetVal/DiracNet,
     exactly as the inlined expression did. *)
  cseDefs = {}; cseDecls = "";
  (* The lc<k>()/dc<k>() net-level CSE runs for BOTH paths: the dressed (collected) lnet builders are
     just as repetitive as the dense σ case (the same Lorentz/projector structures), and the accessors
     only change the GENERATOR's internal sharing — each use copies the shared NetVal/DiracNet, so the
     contracted result (hence the committed kernel) is byte-identical. Previously skipped when dressed,
     which left the dressed generator C++ bloated and its -O0 compile ~2x slower. The SUB-TERM dedup
     (below) stays non-dressed only: a dressed sub-term is keyed by (dn,ln,dch,dsl) not just (dn,ln),
     and the collected path emits ~1 sub-term per net, so there is nothing to merge. *)
  Module[{lTerms = Flatten[lStr], dTerms = Flatten[dStr], lCount, dCount, lMap = <||>, dMap = <||>,
          li = 0, di = 0, lRef, dRef},
    With[{ntT = First@AbsoluteTiming[lCount = Counts[lTerms]; dCount = Counts[dTerms];]},
      ntLog["[prof] CSE Counts (", Length[lTerms], "+", Length[dTerms], " terms): ", ntT, " s"]];
    KeyValueMap[Function[{t, c}, If[c >= 2 && t =!= "NetVal{}" && t =!= "", lMap[t] = "lc" <> bb[li]; li++]], lCount];
    KeyValueMap[Function[{t, c}, If[c >= 2 && t =!= "DiracNet{}" && t =!= "", dMap[t] = "dc" <> bb[di]; di++]], dCount];
    lRef[t_] := If[KeyExistsQ[lMap, t], lMap[t] <> "()", t];
    dRef[t_] := If[KeyExistsQ[dMap, t], dMap[t] <> "()", t];
    cseDefs = Join[
      KeyValueMap[Function[{t, nm}, "const DiracNet& " <> nm <> "(){ static const DiracNet v = " <> t <> "; return v; }"], dMap],
      KeyValueMap[Function[{t, nm}, "const NetVal& " <> nm <> "(){ static const NetVal v = " <> t <> "; return v; }"], lMap]];
    cseDecls = StringJoin[Riffle[Join[
      KeyValueMap[Function[{t, nm}, "const DiracNet& " <> nm <> "();"], dMap],
      KeyValueMap[Function[{t, nm}, "const NetVal& " <> nm <> "();"], lMap]], "\n"]];
    With[{ntT = First@AbsoluteTiming[dStr = Map[dRef, dStr, {2}]; lStr = Map[lRef, lStr, {2}];]},
      ntLog["[prof] CSE ref-rewrite: ", ntT, " s"]];
    If[! hasDressed,
    (* SUB-TERM DEDUP: a net is Σ_b scal_b·contract(dn_b, ln_b); the dense σ^μν projection emits the
       SAME (dn,ln) pair many times with different scalars (measured 48624 sub-terms but only ~22416
       distinct pairs, 2.2x). Merge identical pairs, summing their scalars, and DROP pairs whose scalars
       sum to 0 (σ-commutator cancellations) — so the run-time contraction (the bottleneck) computes
       numeric_value_netval ONCE per distinct pair instead of per occurrence. Output-preserving. *)
    Module[{before = Total[Length /@ dStr], after},
      With[{ntT = First@AbsoluteTiming[
      {dStr, lStr, dscv} = Transpose @ MapThread[Function[{ds, ls, sc},
         Module[{kept = Select[({#[[1, 1]], #[[1, 2]], Total[#[[All, 3]]]} &) /@
                  GatherBy[Transpose[{ds, ls, sc}], #[[{1, 2}]] &], #[[3]] =!= 0 &]},
           If[kept === {}, {{}, {}, {}}, Transpose[kept]]]], {dStr, lStr, dscv}];]},
        ntLog["[prof] sub-term dedup: ", ntT, " s"]];
      after = Total[Length /@ dStr];
      ntLog["[cse] net terms: lnet ", Length[lTerms], "->", Length[lMap], " distinct, dnet ",
            Length[dTerms], "->", Length[dMap], " distinct shared builders; sub-terms ",
            before, "->", after, " after (dn,ln)-pair dedup"]],
      ntLog["[cse] net terms: lnet ", Length[lTerms], "->", Length[lMap], " distinct, dnet ",
            Length[dTerms], "->", Length[dMap], " distinct shared builders (dressed; sub-term dedup skipped)"]]];
  (* each net builder returns the group's LIST of sub-term Dirac chains / Lorentz nets (one function
     per net, summed in main — like the inv body's add(), but type-correct for the numeric backend). *)
  dnetDefs = MapIndexed["std::vector<DiracNet> dnet" <> bb[#2[[1]] - 1] <> "(){ return {" <> StringRiffle[#1, ", "] <> "}; }" &, dStr];
  lnetDefs = MapIndexed["std::vector<NetVal> lnet" <> bb[#2[[1]] - 1] <> "(){ return {" <> StringRiffle[#1, ", "] <> "}; }" &, lStr];
  (* DRESSED: per-net builders for the chain/slot tables, mirroring dnet/lnet — so the giant `dch`/`dsl`
     literals compile on the parallel -O0 units instead of the single -O1 main TU (see unitPre note). *)
  dchDefs = If[hasDressed,
    MapIndexed["std::vector<std::vector<DChainTok>> dch" <> bb[#2[[1]] - 1] <> "(){ return {" <> StringRiffle[#1, ", "] <> "}; }" &, ddch], {}];
  dslDefs = If[hasDressed,
    MapIndexed["std::vector<std::vector<DSlot>> dsl" <> bb[#2[[1]] - 1] <> "(){ return {" <> StringRiffle[#1, ", "] <> "}; }" &, ddsl], {}];
  allDefs = Join[cseDefs, dnetDefs, lnetDefs, dchDefs, dslDefs];
  (* Size-aware unit count (~250 KB/unit, 8..48): the CSE accessors are many but small, so the old
     12-defs/unit rule would emit hundreds of tiny TUs each re-parsing the shared decl header. *)
  nUnits = Min[48, Max[8, Ceiling[StringLength[StringJoin[allDefs]]/250000]]];
  units = Table[unitPre <> StringRiffle[Table[allDefs[[i]], {i, u, Length[allDefs], nUnits}], "\n"] <> "\n",
    {u, 1, Min[nUnits, Length[allDefs]]}];
  (* the shared decl header: net-builder + CSE-accessor forward declarations (the units that call the
     lc<k>()/dc<k>() accessors #include this — emitted ONCE here, not duplicated per unit). *)
  decl = "// GENERATED by MakeNTKernel — do not edit. Numeric net-builder declarations.\n#pragma once\n" <>
    "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n#include <vector>\n" <>
    If[hasDressed, "#include \"numtracer/numeric/numeric_contract.hpp\"\n", ""] <>
    "using namespace numtracer::network;\n" <>
    If[hasDressed, "using namespace numtracer::numeric;\n", ""] <>
    cseDecls <> "\n" <>
    StringJoin[Table["std::vector<DiracNet> dnet" <> bb[i] <> "();\n", {i, 0, nNet - 1}]] <>
    StringJoin[Table["std::vector<NetVal> lnet" <> bb[i] <> "();\n", {i, 0, nNet - 1}]] <>
    If[hasDressed,
      StringJoin[Table["std::vector<std::vector<DChainTok>> dch" <> bb[i] <> "();\n", {i, 0, nNet - 1}]] <>
      StringJoin[Table["std::vector<std::vector<DSlot>> dsl" <> bb[i] <> "();\n", {i, 0, nNet - 1}]], ""];
  (* component-table init: comp[base][mu] = <MPoly builder>, skipping structural zeros. *)
  compInit = StringJoin @ KeyValueMap[Function[{base, comps},
     StringJoin @ MapIndexed[Function[{s, mu},
        If[s === "MPoly(" <> bb[nsym] <> ")", "",
           "  comp[" <> bb[base] <> "][" <> bb[mu[[1]] - 1] <> "] = " <> s <> ";\n"]], comps]], compCpp];
  main = StringJoin[Flatten[{
    "int main(int argc, char** argv){\n",
    "  std::string decor = \"static inline\"; std::string hns = \"" <> nsInner <> "\";\n",
    "  for(int a=1;a<argc;++a){ std::string s=argv[a]; if(s==\"-d\"&&a+1<argc) decor=argv[++a]; else if(s==\"-n\"&&a+1<argc) hns=argv[++a]; }\n",
    "  const int nsym = " <> bb[nsym] <> ";\n",
    "  std::vector<std::array<MPoly,4>> comp(" <> bb[maxBase + 1] <> ", {MPoly(nsym),MPoly(nsym),MPoly(nsym),MPoly(nsym)});\n",
    compInit,
    "  std::vector<std::string> symNames = {" <> StringRiffle[("\"" <> # <> "\"") & /@ symNames, ","] <> "};\n",
    "  std::vector<std::vector<int>> units = {" <>
      StringRiffle[("{" <> StringRiffle[bb /@ #, ","] <> "}") & /@ unitG, ","] <> "};\n",
    "  std::vector<std::vector<DiracNet>> dn = {" <> StringRiffle[Table["dnet" <> bb[i] <> "()", {i, 0, nNet - 1}], ","] <> "};\n",
    "  std::vector<std::vector<NetVal>> ln = {" <> StringRiffle[Table["lnet" <> bb[i] <> "()", {i, 0, nNet - 1}], ","] <> "};\n",
    "  std::vector<std::vector<Cx>> dsc = {" <>
      StringRiffle[Function[scl, "{" <> StringRiffle[Function[s, "Cx{" <> cppNum[Re[s]] <> "," <> cppNum[Im[s]] <> "}"] /@ scl, ","] <> "}"] /@ dscv, ","] <> "};\n",
    (* dressed-numerator chains/slots (symbolic dressing collection): per net, per sub-term. A plain
       sub-term has an empty chain (uses dn[i][j]); a dressed sub-term uses dch[i][j]/dsl[i][j] via
       numeric_value_dressed_netval. Emitted only when the kernel has dressed nets. *)
    If[hasDressed,
      "  std::vector<std::vector<std::vector<DChainTok>>> dch = {" <>
        StringRiffle[Table["dch" <> bb[i] <> "()", {i, 0, nNet - 1}], ","] <> "};\n" <>
      "  std::vector<std::vector<std::vector<DSlot>>> dsl = {" <>
        StringRiffle[Table["dsl" <> bb[i] <> "()", {i, 0, nNet - 1}], ","] <> "};\n", ""],
    "  std::vector<NetVal> lnflat; for(auto &v: ln) for(auto &x: v) lnflat.push_back(x);\n",
    "  auto atomDen = collect_atom_denoms(nsym, lnflat, comp);\n",
    "  for(auto &a: atomDen) a = reduce_units(a, units);  // bare-loop k^2 -> monomial l1^2 -> cancels\n",
    (* NET-PARALLEL contraction: each net's matrix-product contraction is independent (the heavy
       phase). Worker count capped by NT_GEN_MAXW (RAM lever). The serial colour-fold + lowering
       below interns into the shared GlobalEnv (not thread-safe), so only the contraction is threaded.
       Result is order-independent (mp[i] keyed by i), so the kernel is reproducible. *)
    "  std::vector<" <> If[hasDressed, "DPoly", "MPoly"] <> "> mp(" <> bb[nNet] <> ", " <> If[hasDressed, "DPoly(nsym)", "MPoly(nsym)"] <> ");\n",
    "  const bool ntprof = (std::getenv(\"NT_GEN_PROFILE\")!=nullptr);\n",
    "  { unsigned hw=std::thread::hardware_concurrency(); if(!hw)hw=4u;\n",
    "    if(const char* mw=std::getenv(\"NT_GEN_MAXW\")){int v=std::atoi(mw); if(v>0&&(unsigned)v<hw)hw=(unsigned)v;}\n",
    "    std::atomic<int> next{0}; std::mutex lgmx; auto t0=std::chrono::steady_clock::now();\n",
    "    auto work=[&]{ int i; while((i=next.fetch_add(1))<" <> bb[nNet] <> "){\n",
    "        auto a=std::chrono::steady_clock::now(); std::size_t nsub=dn[i].size();\n",
    If[hasDressed,
      "        DPoly acc(nsym);\n" <>
      "        for(std::size_t j=0;j<nsub;++j){\n" <>
      "          DPoly sub = dch[i][j].empty()\n" <>
      "            ? DPoly::fromMPoly(numeric_value_netval(nsym, dn[i][j], ln[i][j], comp, atomDen, units))\n" <>
      "            : numeric_value_dressed_netval(nsym, dch[i][j], dsl[i][j], ln[i][j], comp, atomDen, units);\n" <>
      "          acc = acc + scaleCx(sub, dsc[i][j]); }\n",
      "        MPoly acc(nsym);\n" <>
      "        for(std::size_t j=0;j<nsub;++j) acc = acc + MPoly::constant(nsym, dsc[i][j]) * numeric_value_netval(nsym, dn[i][j], ln[i][j], comp, atomDen, units);\n"],
    "        mp[i]=std::move(acc);\n",
    "        if(ntprof){ double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-a).count();\n",
    "          std::lock_guard<std::mutex> lk(lgmx); std::fprintf(stderr,\"[num] net=%d subterms=%zu mp=%zu  %.1f ms\\n\", i, nsub, mp[i].size(), ms);} } };\n",
    "    int W=(int)std::min<unsigned>(hw,(unsigned)std::max(1," <> bb[nNet] <> "));\n",
    "    std::vector<std::thread> pool; for(int w=1;w<W;++w){ try{pool.emplace_back(work);}catch(const std::system_error&){break;} }\n",
    "    work(); for(auto&t:pool)t.join();\n",
    "    if(ntprof) std::fprintf(stderr,\"[num] contraction phase: %.1f s (W=%u, %d nets)\\n\", std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count(), hw, " <> bb[nNet] <> ");\n",
    "  }\n",
    "  std::vector<SUNNet> colnets = {" <> StringRiffle[colourNets, ","] <> "};\n",
    "  std::vector<numtracer::Cx> colv(" <> bb[nNet] <> "); for(int i=0;i<" <> bb[nNet] <> ";++i) colv[i]=sun_value_cx(colnets[i]);\n",
    "  std::vector<std::vector<int>> groups = {" <>
      StringRiffle[("{" <> StringRiffle[bb /@ #, ","] <> "}") & /@ groups, ","] <> "};\n",
    "  GlobalEnv g;\n",
    "  std::vector<GenProg> progs;\n",
    (* realOnly[gi]: this group's dressing coeff is real, so only Re(trace) is consumed -> emit a
       double trace (no dead imaginary half). Defaults to all-0 (full complex) if not supplied. *)
    "  std::vector<int> realOnly = {" <> StringRiffle[
       If[Length[realOnlyG] === nGrp, (If[TrueQ[#], "1", "0"]) & /@ realOnlyG, Table["0", {nGrp}]], ","] <> "};\n",
    If[hasDressed,
      "  for(size_t gi=0; gi<groups.size(); ++gi){ auto &grp=groups[gi]; DPoly acc(nsym); for(int d: grp) acc = acc + scaleCx(mp[d], colv[d]); progs.push_back(to_genprog(acc, g, realOnly[gi]!=0)); }\n",
      "  for(size_t gi=0; gi<groups.size(); ++gi){ auto &grp=groups[gi]; MPoly acc(nsym); for(int d: grp) acc = acc + mp[d]*MPoly::constant(nsym, colv[d]); progs.push_back(to_genprog(acc, g, realOnly[gi]!=0)); }\n"],
    "  FillFormulas fm;\n",
    "  fm.var = [](int id)->std::string{\n",
    Table["    if(id==" <> bb[i - 1] <> ") return \"" <> varFill[[i]] <> "\";\n", {i, 1, Length[varFill]}],
    "    return \"0.0\"; };\n",
    "  fm.inv = [&](int id)->std::string{ return \"1.0/(\" + mpoly_to_cpp(atomDen[(size_t)id], symNames) + \")\"; };\n",
    (* dressing fill (kind-2 `dress` leaves): the kernel BODY evaluates each dressing atom (where the
       regulators REG::* and the interpolator parameters are in scope) into `dr_<id>` and passes the
       VALUE to fill(); fill just stores it. So fm.dress returns the passed-in argument name. *)
    If[hasDressed,
      "  fm.dress = [](int id)->std::string{ return \"dr_\" + std::to_string(id); };\n", ""],
    "  std::cout << \"// GENERATED by gen_" <> nsInner <> ".cpp — do not edit.\\n\";\n",
    "  std::cout << \"#pragma once\\n#include <cmath>\\n" <> If[complexQ, "#include <complex>\\n", ""] <> "namespace " <> kns <> " { namespace \" << hns << \" {\\n\";\n",
    "  std::cout << \"template<int N> \" << decor << \" double powr(double x){ double r=1.0; for(int i=0;i<N;++i) r*=x; return r; }\\n\";\n",
    "  emit_env_layout(std::cout, g);\n",
    "  std::cout << \"static inline constexpr int nenv = \" << g.syms.size() << \";\\n\";\n",
    "  emit_fill(std::cout, g, \"fill\", \"" <> fillArgSig <> "\", fm, decor);\n",
    (* TRACE-BODY DEDUP: the grouping key is the dressing COEFFICIENT (diagData), finer than the trace
       STRUCTURE — so flows with many Feynman graphs that share a kinematic trace but differ only in
       their dressing/coupling coefficient (e.g. the σ-Yukawa hSigL: 503 groups, only 157 distinct trace
       bodies) emit hundreds of byte-identical trN. Render each trace, key on the name-independent body,
       and emit a duplicate as a one-line forwarder `trN(f){ return trK(f); }` (canonical K = first with
       that body). The body is computed identically once per call site regardless (GCC CSEs the inlined
       identical traces — see ZA4 fusion notes), so this is a pure SOURCE-SIZE/compile win, runtime-neutral.
       Flows with no shared trace structure (ZAqbq1/4/7_147: 108/108 distinct) never hit `seen` ⇒ the
       emitted bytes are unchanged. *)
    "  { std::unordered_map<std::string,std::string> seen; seen.reserve((size_t)" <> bb[nGrp] <> ");\n" <>
    "    for(int i=0;i<" <> bb[nGrp] <> ";++i){\n" <>
    "      std::ostringstream os; emit_cpp(os, progs[i], \"tr\"+std::to_string(i), decor);\n" <>
    "      std::string s = os.str(); std::string body = s.substr(s.find('{'));\n" <>
    "      auto it = seen.find(body);\n" <>
    "      if(it==seen.end()){ seen.emplace(std::move(body), \"tr\"+std::to_string(i)); std::cout << s; }\n" <>
    "      else { const char* rt = (s.find(\"std::complex<double> tr\")!=std::string::npos) ? \" std::complex<double> \" : \" double \";\n" <>
    "        std::cout << decor << rt << \"tr\" << i << \"(const double *f) { return \" << it->second << \"(f); }\\n\"; } } }\n",
    "  std::cout << \"}} // namespace " <> kns <> "::\" << hns << \"\\n\";\n  return 0;\n}\n"}]];
  {pre, units, decl, main}
];


(* ---- whole kernel (boilerplate delegated to FunKit) ------------------------- *)

(* regulator wrappers, prefixed with the user-chosen decorator (default plain "static inline";
   pass e.g. "Decorator" -> "static __host__ __device__ inline" for CUDA-callable kernels) *)
privDefs[decor_] := StringRiffle[
  (decor <> " auto " <> # <> "(const auto &k2, const auto &p2) { return REG::" <> # <> "(k2, p2); }") & /@
    {"RB", "RF", "RBdot", "RFdot", "dq2RB", "dq2RF"}, "\n"];

(* per-trace real/imag accessors for the "RePart" re/im split: a kernel whose VALUE is real but which
   has a COMPLEX trace is assembled as Σ[Re(c)·tr.real() − Im(c)·tr.imag()] — PURE double arithmetic, so
   no complex type survives into the kernel (avoids mixing the trace's std::complex with the support
   complex, and keeps device code real). Overloads pass a real trace straight through (im → 0). *)
ntReImDefs[decor_] := StringRiffle[{
  decor <> " double ntRe(double x) { return x; }",
  "template <class T> " <> decor <> " double ntRe(const T &z) { return z.real(); }",
  decor <> " double ntIm(double) { return 0.0; }",
  "template <class T> " <> decor <> " double ntIm(const T &z) { return z.imag(); }"}, "\n"];

(* ---- emit helpers: support `using`s, namespace wrapping, runtime include -------------------
   A generated kernel pulls its math helpers (complex, powr, pow, sqrt, fma) from a configurable
   SUPPORT namespace and is wrapped in a configurable HOST namespace. The defaults (numtracer /
   numtracer/codegen/runtime.hpp) make the emitted code self-contained against NumTracer's own
   headers; a consumer that already provides equivalents (e.g. DiFfRG) points the codegen at
   them via "SupportNamespace" / "KernelNamespace" / "RuntimeInclude". *)
ntSupportUsings[sns_] := DeleteDuplicates[{
  "using namespace " <> sns <> ";", "using namespace " <> sns <> "::compute;",
  "using namespace numtracer;"}];
(* the loop-independent `constant(p, k, dressings...)` function. DiFfRG flat-adds its return to the
   integral (constExpr second arg of MakeKernel); NumTracer emits the same. A 0 constant needs no
   namespace usings (and stays byte-identical to the pre-Constant emission); a nonzero expr may call
   compute helpers (powr/pow/…) or the support API, so it gets the same usings the kernel body does. *)
ntConstFn[constExpr_, decor_, constParams_, sns_] := FunKit`MakeCppFunction[constExpr,
  "Name" -> "constant", "Prefix" -> decor, "Return" -> "auto", "CodeParser" -> "Cpp",
  "Parameters" -> constParams,
  "Body" -> If[MatchQ[constExpr, 0 | 0.], "", StringRiffle[ntSupportUsings[sns], "\n"]]];
ntWrapBody[kns_, classStr_, name_] := If[kns === None || kns === "",
  {classStr},
  {"namespace " <> kns <> "\n{", classStr, "}", "using " <> kns <> "::" <> name <> ";"}];
ntRuntimeIncludes[runInc_] := If[runInc === None || runInc === "", {}, {runInc}];
(* the dressing-parameter type: Automatic -> `const auto&` (fully generic, self-contained);
   else the given concrete type string (e.g. a consumer's interpolator type). *)
ntDressType[dressTy_] := If[dressTy === Automatic, "auto", dressTy];

(* ---- direct (single-header) backend, the MakeNTKernel "Dense" path ---------------
   "Dense" (brute-force numeric DTensor oracle): emits `numtracer::dense::` builder calls read via
   `.scalar_value().re`, used to cross-check the numeric matrix-product backend. *)
Options[mkDirectKernel] = {"Name" -> "nt_kernel",
  "Dressings" -> {}, "ScalarParams" -> {}, "Backend" -> "Dense", "Decorator" -> "static inline", "AngleDefs" -> {},
  "RuntimeInclude" -> "numtracer/codegen/runtime.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "numtracer_kernels", "SupportNamespace" -> "numtracer",
  "DressingType" -> Automatic, "RegulatorAlias" -> False, "Constant" -> 0.};
(* "Constant" -> expr (default 0.): loop-independent flat-added term -> body of constant(); see
   the matching option on mkGenerateKernel / MakeNTKernel. *)
(* "Decorator" -> "<prefix>": the function prefix emitted on kernel/constant and the regulator
   wrappers (e.g. "static __host__ __device__ inline" so the kernel is CUDA-device-callable).
   The default keeps the emitted header byte-identical to before the option existed. *)

(* Single-header kernel emitter for the ET / Dense backends (NOT the build-time generator path —
   that is mkGenerateKernel). Algorithm:
     1. unpack the NTKernel `k0` and resolve options (name, dressings, backend, namespaces, decorator).
     2. rebuild the fundamental-symbol env (`env`/`mask`) by packing each diagram's momentum components,
        and build the per-symbol fill formulas.
     3. walk the diagrams (compileDiagram / compileDiagramDense) to collect the kernel-body preamble
        (the per-component trace declarations) and the integrand `Σ coeff_i × trace_i`.
     4. hand the preamble + integrand + fill to FunKit's COEN (MakeCppClass) to emit the C++ class with
        kernel()/constant() methods, and write the header (write-if-changed). *)
mkDirectKernel[NTKernel[k0_], file_, OptionsPattern[]] := Module[
  {name, dress, scalarParams, backend, complexQ, args, frame, env, nenv, nc, mask, fills, preamble, integrand,
   kernelParams, constParams, mkParam, kernelFn, constFn, classStr, header, decor, angleDefs, angleDecls,
   kns, sns, runInc, extraInc, interpTy, k},
  Needs["FunKit`"];
  (* The DENSE backend cannot fold an eager dressed sum (it contracts numerically, no runtime-dressing
     atoms), so if the kernel was built with symbolic dressing collection (ntDressedNum tokens), expand
     them back to the distributed Dirac structure sum and RE-DISTRIBUTE each diagram (collection OFF) —
     reproducing exactly what the dense oracle saw before collection. The numeric backend keeps the
     collected form; only this dense path re-distributes. *)
  k = If[FreeQ[k0["Diagrams"], _ntDressedNum], k0,
    <|k0, "Diagrams" -> Flatten[redistDiagram /@ k0["Diagrams"]]|>];
  name  = OptionValue["Name"];  dress = OptionValue["Dressings"];
  scalarParams = OptionValue["ScalarParams"];  (* loop-independent scalar doubles threaded into the signature *)
  angleDefs = OptionValue["AngleDefs"];  (* symbolic kinematic angles (cosl1pj) defined once in the body *)
  backend = OptionValue["Backend"];  decor = OptionValue["Decorator"];
  kns = OptionValue["KernelNamespace"];  sns = OptionValue["SupportNamespace"];
  runInc = OptionValue["RuntimeInclude"];  extraInc = OptionValue["ExtraIncludes"];
  interpTy = ntDressType[OptionValue["DressingType"]];
  (* a Dirac-vertex projector carries an `i`: the kernel computes complex intermediates whose
     SUM is real (FORM shows this). In the dense backend we keep the momentum trace tokens
     complex and take Re of the assembled integrand. (The et backend is dropped for such
     vertices — its symbolic type is too big to compile.) *)
  complexQ = backend === "Dense" && ! FreeQ[k["Diagrams"], Complex];
  args  = k["Args"];  frame = k["Frame"];  env = k["Env"];  nenv = k["NEnv"];
  (* Vestigial compile-chain salt: every SU(N) head now carries its own rank N, so the dense/inv
     builders read it from the head and never consult this positional arg (kept only to avoid
     re-threading the recursion). *)
  nc = 0;

  mask = Association @ KeyValueMap[#1 -> frameMask[resolveComponents[#1, frame]] &, env];

  (* Pack renv: each momentum keeps only its structurally-nonzero components (popcount[mask]
     slots instead of a fixed 4); the 1/q^2 slots follow all component slots. This repack is
     LOCAL to the Dense path — buildEnv (shared with the numeric backend, which never emits
     renv) keeps its 4-per-Base layout. The dense builders read the k-th present component (in
     Lorentz order) from renv[Base+k]; `mask` still carries the Lorentz placement. *)
  Module[{base = 0, inv},
    env = Association @ Map[
      Function[q, q -> <|env[q], "Base" -> With[{b = base}, base += DigitCount[mask[q], 2, 1]; b]|>],
      Keys[env]];
    inv = base;
    env = Association @ Map[
      Function[q, q -> <|env[q], "Inv" -> If[env[q]["Inv"] === None, None, inv++]|>],
      Keys[env]];
    env = Association @ Map[                                  (* spatial 1/|q⃗|² slots (finite-T E/M) *)
      Function[q, q -> <|env[q], "InvS" -> If[Lookup[env[q], "InvS", None] === None, None, inv++]|>],
      Keys[env]];
    nenv = inv];

  (* renv fill: each structurally-nonzero momentum component (packed), then the 1/q^2 slots *)
  fills = Flatten @ {
    KeyValueMap[Function[{q, info}, Module[{comps = resolveComponents[q, frame], b = info["Base"], kk = 0},
      Table[If[comps[[i + 1]] =!= 0 && comps[[i + 1]] =!= 0.,
        "renv[" <> ToString[b + kk++] <> "] = " <> FunKit`CppForm[comps[[i + 1]]] <> ";", Nothing], {i, 0, 3}]]], env],
    KeyValueMap[Function[{q, info}, If[info["Inv"] =!= None,
      "renv[" <> ToString[info["Inv"]] <> "] = 1.0 / (" <> FunKit`CppForm[normSqExpr[q, frame]] <> ");", Nothing]], env],
    KeyValueMap[Function[{q, info}, If[Lookup[info, "InvS", None] =!= None,
      "renv[" <> ToString[info["InvS"]] <> "] = 1.0 / (" <> FunKit`CppForm[spatialNormSqExpr[q, frame]] <> ");", Nothing]], env]};

  (* walk the diagrams: collect preamble lines + build the integrand with placeholders.
     Each diagram is ONE eager-summed contraction per component (no monomial blow-up). *)
  preamble = {}; integrand = 0;
  MapIndexed[Function[{diag, di}, Module[{decls, coeff, toks},
    {decls, coeff, toks} = compileDiagramDense[diag, ToString[di[[1]] - 1], frame, env, mask, nc, complexQ];
    preamble = Join[preamble, decls];
    integrand += coeff * Times @@ toks]], k["Diagrams"]];
  (* complex dense kernel: the integrand is left COMPLEX (a Dirac-vertex projector carries an
     `i`; the physical flow is its real part). We do NOT take Re in Mathematica — wrapping the
     CSE'd integrand in Re[] triggers a pathological ComplexExpand/Simplify. Instead the kernel
     returns std::complex<double> and the caller takes std::real(...) (which also passes a plain
     double through unchanged, so real kernels are unaffected). *)

  (* kinematic angle defs (kept symbolic in the dressing): emit once as named temporaries, exactly
     like the Generate path, so a Dense baseline for a symmetric-point flow (cosl1pj symbols) compiles. *)
  angleDecls = ("const double " <> SymbolName[First[#]] <> " = " <> cppFlat[Last[#]] <> ";") & /@ angleDefs;
  preamble = StringRiffle[Join[
    ntSupportUsings[sns],
    {If[complexQ, "using std::complex;", Nothing],
     "double renv[" <> ToString[Max[nenv, 1]] <> "] = {};"}, angleDecls, fills, preamble], "\n"];

  mkParam[nm_, ty_] := <|"Name" -> If[StringQ[nm], nm, SymbolName[nm]], "Type" -> ty, "Const" -> True, "Reference" -> True|>;
  (* scalar doubles (etaPiL, d1V, rhoL, ...) sit between args and dressings, matching the order
     DiFfRG's integrator forwards them (k, scalars, interpolators); they appear in the kernel /
     constant BODY (propagator masses, anomalous-dim constants) but never in the fill. *)
  kernelParams = Join[mkParam[#, "double"] & /@ args, mkParam[#, "double"] & /@ scalarParams, mkParam[#, interpTy] & /@ dress];
  constParams  = Join[mkParam[#, "double"] & /@ Select[args, # === Global`p || # === Global`k &],
                      mkParam[#, "double"] & /@ scalarParams,
                      mkParam[#, interpTy] & /@ dress];

  kernelFn = FunKit`MakeCppFunction[integrand, "Name" -> "kernel", "Prefix" -> decor,
    "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> kernelParams, "Body" -> preamble];
  constFn = ntConstFn[OptionValue["Constant"], decor, constParams, sns];

  classStr = FunKit`MakeCppClass["TemplateTypes" -> {"REG"}, "Name" -> name,
    "MembersPublic" -> If[TrueQ[OptionValue["RegulatorAlias"]],
      Prepend[{kernelFn, constFn}, "using Regulator = REG;"], {kernelFn, constFn}],
    "MembersPrivate" -> {privDefs[decor]}];

  header = FunKit`MakeCppHeader[
    (* the dense umbrella pulls in the axis planner + dirac_data + sun_data + cx transitively *)
    "Includes" -> Join[extraInc, ntRuntimeIncludes[runInc], {"numtracer/dense/dtensor.hpp"}],
    "Body" -> ntWrapBody[kns, classStr, name]];

  (* header is already clang-formatted by MakeCppHeader; write-if-changed *)
  If[FileExistsQ[file] && Import[file, "Text"] === header,
    Print["unchanged: ", file],
    Export[file, header, "Text"]; Print["wrote: ", file]];
  file
];

(* ---- numeric matrix-product kernel, the MakeNTKernel "Numeric" path: emit the build-time
        generator program + the kernel header that calls its output. ---------------------- *)

Options[mkGenerateKernel] = {"Name" -> "nt_inv_kernel", "Namespace" -> Automatic,
  "Dressings" -> {}, "ScalarParams" -> {}, "IncludeDir" -> Automatic,
  "RunGenerator" -> True, "FullParallel" -> False, "AngleDefs" -> {}, "CrossTraceCSE" -> False,
  "GlobalCollect" -> True,  (* default ON: groups diagrams by dressing coeff + folds colour
    numerically so COEN collects the factored coefficient across diagrams (FORM-style); ~30% runtime
    cut on the dense quark-loop vertices, generation cost unchanged. Pass False to opt out per flow. *)
  "NumericContract" -> False, "Components" -> Automatic, "SymbolDefs" -> <||>,
  "Decorator" -> "static inline",
  "RuntimeInclude" -> "numtracer/codegen/runtime.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "numtracer_kernels", "SupportNamespace" -> "numtracer",
  "DressingType" -> Automatic, "RegulatorAlias" -> False, "RealProbe" -> True,
  "PruneRealTraces" -> False, "Constant" -> 0.};
(* The dressing-collection path is driven by the `ntDressedNum` tokens NumTrace emits under
   "DressingCollection" -> True (no MakeNTKernel option needed): mkGenerateKernel auto-detects them and
   routes to the DPoly generator branch. *)
(* "Constant" -> expr (default 0.): the loop-INDEPENDENT piece that DiFfRG flat-adds to the
   integral — emitted verbatim as the body of `constant(p, k, dressings...)`, exactly like the
   constExpr second argument of DiFfRG's MakeKernel. A plain Mathematica expression in p/k and the
   dressing names (e.g. ZA[p] -> ZA(p)); NOT an NTKernel. Left at 0. the constant returns 0.. *)
(* "PruneRealTraces" -> False (default OFF): emit a `double` trace for any diagram group whose dressing
   coefficient is real (only Re(trace) is consumed), skipping the dead imaginary half. CORRECT for the
   kernel in every verdict, but HAZARDOUS with the probe: the probe runs on the generated traces and
   verifies Im(integrand)≈0; dropping a real-coeff group's imaginary half removes a term the probe
   needs, which can leave a non-cancelling residual and misclassify a real flow as Complex (blocking
   the double-kernel emission). Safe to enable only when the probe is off OR after per-flow validation,
   or once a post-probe trace regeneration is wired. Left off pending measurement (GCC often already
   DCEs the dead half of small inlined traces). *)
(* "RealProbe" -> True (default): when the syntactic complexQ trips (some diagram coeff carries an `i`),
   compile+run a probe over the JUST-generated real traces to test whether Im(integrand) actually
   vanishes (projector-i × colour-i often cancel to a real value that Mathematica can't see through the
   opaque trace symbols). If it vanishes, re-emit a REAL (double) kernel — losslessly. Set False to skip
   the probe and always keep the complex+consumer-Re path. *)
(* "NumericContract" -> True: the numeric (matrix-product) backend (task #22) — contract each diagram
   numerically via et/numeric (γ/metric/projector numeric, momenta symbolic), no sp-invariant
   reduce/rebase/ibp. "Components" -> <|mom -> {e0,e1,e2,e3}, ...|> gives each momentum's 4 components
   as expressions (partially numeric / partially symbolic); Automatic falls back to the kinematic
   frame polynomialised. "SymbolDefs" -> <|sym -> expr|> gives the C++ fill for any DERIVED symbol
   (e.g. sin1 -> Sqrt[1-cos1^2]); plain free symbols are taken to be kernel arguments. *)
(* Standalone-output options (defaults make the emitted code self-contained against NumTracer's
   own headers, with no mention of any downstream consumer). A consumer that supplies its own
   support API points the codegen at it via these:
     "RuntimeInclude" -> "<hdr>" | None : the support header #included first, providing `complex`
        and `compute::{powr,pow,sqrt,fma}` (default numtracer/codegen/runtime.hpp; None to omit).
     "ExtraIncludes" -> {"a.hpp", ...}  : extra #includes prepended ahead of everything.
     "KernelNamespace" -> "ns" | None   : namespace wrapping the kernel class AND the generated
        trace functions (default "numtracer_kernels"; None emits at the includer's scope).
     "SupportNamespace" -> "ns"         : where `complex`/`compute` are looked up via `using`
        (default "numtracer").
     "DressingType" -> Automatic | "T"  : dressing-parameter type; Automatic emits `const auto&`
        (fully generic), or give a concrete type string. *)
(* "RunGenerator" -> False: emit the generator sources only, without compiling/running them
   (the committed traces header is left untouched). *)
(* "Decorator" -> "<prefix>": the function prefix on EVERY emitted function — kernel/constant,
   the regulator wrappers, and (via the generator's runtime `-d` flag) fill/trN/powr in the
   straight-line header — e.g. "static __host__ __device__ inline" makes the whole kernel
   CUDA-device-callable. Default keeps the emitted bytes identical (kernel md5 invariant). *)
(* "CrossTraceCSE" -> True: lower all diagram trace polynomials through ONE shared CSE program
   (emitted as trace_all(f, t[]); the kernel fills t[] once and reads t[d]) so subexpressions are
   shared ACROSS traces — a fused kernel. Pays off when many traces share intermediates; default
   False keeps one independent trN() per trace. *)
(* "AngleDefs" -> {sym -> expr, ...}: kinematic angle symbols the dressing keeps SYMBOLIC
   (e.g. cosl1p2 -> (-cos1 + Sqrt[3-3 cos1^2] cos2)/2). Emitted ONCE as `const double sym = ...;`
   in the kernel body so a shared sub-expression (the sqrt) is computed once rather than inlined
   per occurrence. *)
(* "FullParallel" -> True passes `-p` to the generator: heavy nets are reduced/rebased CONCURRENTLY
   (faster codegen, higher peak RAM) instead of one-at-a-time. The emitted kernel is identical. *)

mkGenerateKernel::genfail = "Generator compile/run failed: `1`";
mkGenerateKernel::crosscseComplex = "CrossTraceCSE is incompatible with a complex flow: trace_all \
fills a real `double tarr[]`, which truncates any complex trace's imaginary part, so the RePart \
re/im split would silently read ntIm(tarr[i])==0. Disabling CrossTraceCSE for this flow \
(GlobalCollect subsumes it).";
mkGenerateKernel::emptynets = "Flow `1` produced no generator nets (nets=`2`, groups=`3`) — nothing \
to emit. Aborting instead of writing a placeholder kernel. (A diagram-build guard such as \
disconnectmix may have dropped everything, or NumTrace returned no usable diagrams.)";
mkGenerateKernel::scalarleak = "Diagram `1`: a non-numeric factor `2` reached the generator scalar \
coefficient (a Lorentz tensor that was not resolved by the net builder, e.g. an un-anchored metric \
contraction). It would be emitted as undeclared C++. Aborting; fix the net build (compileTInv) so \
the contraction folds numerically.";

(* Locate the NumTracer C++ headers for the generator compile. In-tree the package sits in
   <repo>/numtracer/mathematica/ with the headers in ../include; installed (e.g. under
   $UserBaseDirectory/Applications/NumTracer) the header location is recorded at configure time
   in the CMake-generated sibling NumTracerPaths.m. Overridable per call ("IncludeDir" option)
   or globally (NUMTRACER_INCLUDE_DIR environment variable). *)
resolveIncludeDir::nodir = "Cannot locate the NumTracer C++ headers (tried the NUMTRACER_INCLUDE_DIR environment variable, the in-tree ../include, the installed NumTracerPaths.m record, and ~/.local/share/NumTracer/include). Pass \"IncludeDir\" -> dir to MakeNTKernel.";
resolveIncludeDir[] := Module[{envDir, dir, pathsFile},
  envDir = Environment["NUMTRACER_INCLUDE_DIR"];
  If[StringQ[envDir] && DirectoryQ[envDir], Return[envDir]];
  dir = FileNameJoin[{DirectoryName[$NumTracerDirectory], "include"}];
  If[DirectoryQ[dir], Return[dir]];
  pathsFile = FileNameJoin[{$NumTracerDirectory, "NumTracerPaths.m"}];
  If[FileExistsQ[pathsFile],
    Get[pathsFile];
    If[StringQ[$NumTracerInstalledIncludeDir] && DirectoryQ[$NumTracerInstalledIncludeDir],
      Return[$NumTracerInstalledIncludeDir]]];
  dir = FileNameJoin[{$HomeDirectory, ".local", "share", "NumTracer", "include"}];
  If[DirectoryQ[dir], Return[dir]];
  Message[resolveIncludeDir::nodir]; Abort[]];

(* ---- semantic complexQ: probe whether the imaginary part actually vanishes -------------------
   The syntactic `complexQ = !FreeQ[Diagrams, Complex]` only sees that SOME diagram coefficient carries
   an `i` (projector i, non-abelian/quark colour f^abc T^bT^c = (iN/2)T^a). Those i's frequently pair up
   (i·i = -1) across diagrams so the assembled flow is exactly REAL — but the cancellation involves the
   concrete trace VALUES, and at the Mathematica stage the traces are opaque generated C++ symbols, so
   it can't be seen there. The just-generated C++ traces CAN see it: this compiles+runs a tiny probe that
   evaluates Im(integrand) over random frames with smooth real stub dressings/regulators. Returns True iff
   the imaginary part vanishes (|Im| <= tol·|Re|) over all sampled points — in which case the caller
   re-emits a real (double) kernel losslessly. Conservative: any failure (compile/run/no-points) -> False
   (keep the complex kernel). *)
Options[numericImagProbeRealQ] = {"NPoints" -> 4000, "Tol" -> 1.*^-9};
numericImagProbeRealQ[integrand_, args_, fillArgs_, angleDefs_, angleDecls_, nsHome_, headerFile_,
    drTable_:<||>, opts : OptionsPattern[]] := Module[
  {keepHeads, keepSyms, seedOf, argComb, stub, probeFull, probeProj, probeParams, probePre, fnFull, fnProj,
   drDecls, drFillArgs, randDecls, callArgs, src, cppFile, bin, rc, out, np, tol, tracesDir, cxx, parsed, distOf},
  np = OptionValue["NPoints"];  tol = OptionValue["Tol"];
  tracesDir = DirectoryName[headerFile];  cxx = resolveGenCxx[];
  (* GENERAL stubbing: replace EVERY external real-valued atom — any dressing (any arity), any
     regulator/support function, any named constant — with `ntStub(seed_head, hash(args))`, an
     INDEPENDENTLY-SEEDED pseudo-random real per head. "External" = head is not a structural math
     operation, not a frame arg, not a known constant, and not a trace token (those are raw C++ strings,
     left untouched). Independence is essential: a single shared stub would make `dress1[x]-dress2[x]`
     collapse to 0 and mask a real surviving imaginary part. Distinct heads -> distinct seeds; distinct
     arguments -> distinct hashes -> distinct values, so any genuine imaginary part is exposed for
     ARBITRARY dressings. *)
  keepHeads = Alternatives[Plus, Times, Power, Rational, Sqrt, Sin, Cos, Tan, Cot, Sec, Csc, Exp, Log,
     Abs, Sign, ArcTan, ArcSin, ArcCos, Sinh, Cosh, Tanh, Max, Min, Floor, Ceiling, Mod, Complex, List,
     Global`ntStub];
  keepSyms = Join[args, First /@ angleDefs, {Pi, E, EulerGamma, Degree, GoldenRatio}];
  seedOf[h_] := N[Mod[Hash[SymbolName[h]], 100003] + 7];      (* distinct seed per external head *)
  argComb[a_List] := Total[MapIndexed[#1 * N[GoldenRatio]^(#2[[1]] - 1) &, a]];  (* 0 for a constant *)
  stub[e_] := Module[{x},
     x = e //. (h_Symbol)[a___] /; FreeQ[keepHeads, h] && ! MemberQ[args, h] :>
        Global`ntStub[seedOf[h], argComb[{a}]];
     x /. (s_Symbol) /; ! MemberQ[keepSyms, s] && FreeQ[keepHeads, s] &&
        Context[s] =!= "System`" && s =!= Global`ntStub :> Global`ntStub[seedOf[s], 0.]];
  probeFull = stub[integrand];
  probeProj = probeFull /. Complex[a_, b_] :> a;        (* the candidate real projection *)
  probeParams = (<|"Name" -> SymbolName[#], "Type" -> "double", "Const" -> True, "Reference" -> True|>) & /@ args;
  (* DRESSED kernels: the generated `fill()` takes one extra `double dr_<id>` per dressing atom (the
     kernel body computes each atom from regulators/interpolators and passes the VALUE). The probe has
     none of that runtime in scope, so it must compute each atom with the SAME pseudo-random stubbing
     used for the integrand and pass the values — otherwise the probe `.cpp` fails to compile (too few
     args to fill) and the kernel is conservatively kept COMPLEX, silently losing the lossless RePart
     double-kernel emission (the imaginary-half DCE that makes the distributed baseline fast). Atoms can
     reference the derived angles, so the angle decls precede them. *)
  drDecls = KeyValueMap[Function[{id, atom},
     "const double dr_" <> ToString[id] <> " = " <> cppFlat[stub[atom]] <> ";"], drTable];
  drFillArgs = ("dr_" <> ToString[#]) & /@ Sort[Keys[drTable]];
  probePre = StringRiffle[Join[
    angleDecls,
    {"double fenv[(" <> nsHome <> "::nenv) > 0 ? (" <> nsHome <> "::nenv) : 1];"},
    drDecls,
    {nsHome <> "::fill(fenv, " <> StringRiffle[Join[SymbolName /@ fillArgs, drFillArgs], ", "] <> ");"}], "\n"];
  fnFull = FunKit`MakeCppFunction[probeFull, "Name" -> "probe_full", "Prefix" -> "static inline",
    "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> probeParams, "Body" -> probePre];
  fnProj = FunKit`MakeCppFunction[probeProj, "Name" -> "probe_proj", "Prefix" -> "static inline",
    "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> probeParams, "Body" -> probePre];
  distOf[a_] := Which[StringContainsQ[SymbolName[a], "cos"], "Uc",
                      StringContainsQ[SymbolName[a], "phi"], "Uph", True, "U"];
  randDecls = StringRiffle[("double " <> SymbolName[#] <> " = " <> distOf[#] <> "(rng);") & /@ args, " "];
  callArgs = StringRiffle[SymbolName /@ args, ", "];
  src = StringJoin[
    "#define __host__\n#define __device__\n",
    "#include <complex>\n#include <cmath>\n#include <random>\n#include <cstdio>\n",
    "#include \"" <> FileNameTake[headerFile] <> "\"\n",
    "template<int N> static inline double powr(double x){ double r=1.0; int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?1.0/r:r; }\n",
    "using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;\n",
    "static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}\n",
    "template<class T> using complex = std::complex<T>;\n",
    (* independently-seeded pseudo-random real in [0.4,0.9): same (seed,arg) -> same value (a dressing is
       a function), distinct (seed,arg) -> independent value, so no two dressings or arguments collide. *)
    "static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }\n",
    fnFull, "\n", fnProj, "\n",
    "int main(){ std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);\n",
    "  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0; long ok=0;\n",
    "  for(int n=0;n<" <> ToString[np] <> ";++n){ " <> randDecls <> "\n",
    "    std::complex<double> f = probe_full(" <> callArgs <> "); std::complex<double> pj = probe_proj(" <> callArgs <> ");\n",
    "    double im=std::imag(f), re=std::real(f), df=std::abs(f-pj);\n",
    (* PER-POINT relative measures (mrim, mrdiff): a global max|Im|/max|Re| can let a localized
       imaginary part hide behind a large |Re| at some OTHER point (catastrophic cancellation). The
       +1 floor degrades gracefully to an absolute test when |Re|/|f| are small. The verdict keys on
       these; the absolute trio is kept only for the log. *)
    "    if(std::isfinite(im)&&std::isfinite(re)&&std::isfinite(df)){ mim=std::max(mim,std::fabs(im)); mdiff=std::max(mdiff,df); mre=std::max(mre,std::fabs(re));\n",
    "      mrim=std::max(mrim, std::fabs(im)/(std::fabs(re)+1.0)); mrdiff=std::max(mrdiff, df/(std::abs(f)+1.0)); ++ok; } }\n",
    "  std::printf(\"%.10e %.10e %.10e %.10e %.10e %ld\\n\", mim, mdiff, mre, mrim, mrdiff, ok); return 0; }\n"];
  cppFile = FileNameJoin[{$TemporaryDirectory, "ntprobe_" <> StringReplace[nsHome, {":" -> "_"}] <> ".cpp"}];
  bin = StringReplace[cppFile, ".cpp" -> ""];
  Export[cppFile, src, "Text"];
  rc = Run[cxx <> " -std=c++20 -O1 -w -I '" <> tracesDir <> "' '" <> cppFile <> "' -o '" <> bin <> "' 2> '" <> bin <> ".cerr'"];
  If[rc =!= 0, ntLog["[probe] compile failed (rc=", rc, ") — keeping complex kernel (conservative)"]; Return["Complex"]];
  Run["'" <> bin <> "' > '" <> bin <> ".out'"];
  out = If[FileExistsQ[bin <> ".out"], Import[bin <> ".out", "Text"], ""];
  parsed = Quiet@Check[
     ToExpression[StringReplace[#, {"e+" -> "*^", "e-" -> "*^-", "e" -> "*^"}]] & /@ StringSplit[StringTrim[out]],
     $Failed];
  If[! MatchQ[parsed, {_?NumericQ, _?NumericQ, _?NumericQ, _?NumericQ, _?NumericQ, _?NumericQ}] || parsed[[6]] < 1,
    ntLog["[probe] run produced no usable points — keeping complex kernel"]; Return["Complex"]];
  ntLog["[probe] over ", Round[parsed[[6]]], " pts:  max|Im|=", ScientificForm[parsed[[1]], 3],
        "  max|full-proj|=", ScientificForm[parsed[[2]], 3], "  max|Re|=", ScientificForm[parsed[[3]], 3],
        "  rel|Im|=", ScientificForm[parsed[[4]], 3], "  rel|full-proj|=", ScientificForm[parsed[[5]], 3]];
  (* Three-way verdict from the C++ evaluation (which resolves every complex multiplication):
       "Complex"  Im survives                          -> genuinely complex, keep it.
       "Pure"     Im=0 AND Complex->Re projection exact -> drop imaginary coeffs (clean real arithmetic).
       "RePart"   Im=0 but projection differs           -> a trace is itself complex; the value is real
                                                           but only `.real()` of the full complex result
                                                           is correct, so wrap the return in a C++ real part. *)
  (* verdict on the PER-POINT relative measures (mrim, mrdiff): no global-scale inflation. *)
  Which[parsed[[4]] > tol, "Complex",
        parsed[[5]] <= tol, "Pure",
        True, "RePart"]];

(* ---- group-diagonal dressing fold: SUNPoly via the validated C++ engine ---------------------
   Each diag-dressed colour-net STRING (carrying SUN::diag{Fund,Adj}(...,dr) factors) is folded by
   sun_value_dressed in a tiny build-time program (the same emit/compile/run seam as the imaginary
   probe), returning per net a list of terms {coeffRe, coeffIm, {{dr, component}...}}. Reuses the
   numeric engine verbatim, so the per-component colour weights are byte-identical to the typed-out
   SU(N) tables — no Mathematica reimplementation of the group algebra. *)
diagColPolys[colnetStrs_, includeDir_] := Module[
  {cxx = resolveGenCxx[], src, cppFile, bin, rc, out, lines, res = {}, cur = Null, num},
  num[s_] := ToExpression[StringReplace[s, {"e+" -> "*^", "e-" -> "*^-", "e" -> "*^"}]];
  src = StringJoin[
    "#include \"numtracer/network/sun_net.hpp\"\n#include <cstdio>\n#include <vector>\n",
    "using namespace numtracer; using namespace numtracer::network;\n",
    "int main(){\n  std::vector<SUNNet> nets = {" <> StringRiffle[colnetStrs, ", "] <> "};\n",
    "  for(std::size_t n=0;n<nets.size();++n){\n",
    "    SUNPoly p = sun_value_dressed(nets[n]);\n",
    "    std::printf(\"NET %zu %zu\\n\", n, p.size());\n",
    "    for(const auto& t : p){\n",
    "      std::printf(\"T %.17g %.17g %zu\", t.coeff.re, t.coeff.im, t.dress.size());\n",
    "      for(const auto& d : t.dress) std::printf(\" %d %d\", d.first, d.second);\n",
    "      std::printf(\"\\n\"); } }\n  return 0;\n}\n"];
  cppFile = FileNameJoin[{$TemporaryDirectory, "ntdiagpoly.cpp"}];
  bin = StringReplace[cppFile, ".cpp" -> ""];
  Export[cppFile, src, "Text"];
  rc = Run[cxx <> " -std=c++20 -O1 -w -I '" <> includeDir <> "' '" <> cppFile <> "' -o '" <> bin <>
           "' 2> '" <> bin <> ".cerr'"];
  If[rc =!= 0,
    Print["[diagpoly] compile failed (rc=", rc, "):\n", Quiet@Check[Import[bin <> ".cerr", "Text"], ""]];
    Abort[]];
  Run["'" <> bin <> "' > '" <> bin <> ".out'"];
  out = If[FileExistsQ[bin <> ".out"], Import[bin <> ".out", "Text"], ""];
  lines = Select[StringSplit[StringTrim[out], "\n"], # =!= "" &];
  Do[Module[{tk = StringSplit[ln]},
     Which[
       tk[[1]] === "NET", If[cur =!= Null, AppendTo[res, cur]]; cur = {},
       tk[[1]] === "T",
         Module[{re = num[tk[[2]]], im = num[tk[[3]]], m = ToExpression[tk[[4]]]},
           AppendTo[cur, {re, im, If[m === 0, {}, Partition[ToExpression /@ tk[[5 ;; 4 + 2 m]], 2]]}]]]],
    {ln, lines}];
  If[cur =!= Null, AppendTo[res, cur]];
  res];

(* The invariant-basis GENERATION path (the default numeric backend). Algorithm:
     1. reset the per-generation registries ($ctCache, resetDiagDr, resetDr) and unpack the NTKernel.
     2. per diagram, build the Lorentz/colour nets: splitColourGroupsInv groups branches by colour
        structure; a plain diagram lowers to a LorentzNet, a dressed-numerator one to a DPoly chain
        (the two paths the generator contracts differently). CSE-share repeated net/trace builders.
     3. emit a build-time C++ generator program (emitNumericGenerator) — a main TU plus per-net units —
        that runs reduce -> rebase -> lower at codegen time and PRINTS the straight-line traces header.
     4. compile + run that generator (resolveGenCxx), probe whether the result is real/complex to pick
        the kernel's number type, then emit and write (write-if-changed) the kernel header that fills
        the fundamental symbols and calls the generated trN(f). *)
mkGenerateKernel[NTKernel[k_], genFile_, kernelFile_, headerFile_, OptionsPattern[]] := Module[
  {name, ns, dress, scalarParams, args, frame, env, nc, mask, ncomp, fillArgs, fillArgSig, invNets, invRest, g,
   colourNets, gcol, preamble, integrand, kernelParams, constParams, mkParam, kernelFn, constFn,
   classStr, header, hdrInc, incDir, genPre, genUnits, genDecl, genMain, declFile, unitFiles, genSrc, bin, run,
   hasFund, complexQ, colDecls, colToks, angleDefs, angleDecls, crossCSE, traceRef, nGrp, decor,
   kns, sns, runInc, extraInc, interpTy, nsHome, gc, symDefs = <||>, dmono = {}, atomStrs = {},
   groupCombos = {}, groupContribs = {}, realOnlyG = {}, dressedIdx = {}, diagTokExpr = {}, diagDressDims = <||>,
   factorNets = {}, lorFacOf = {}, pGroupOf = <||>, nAdd = 0, factorCompOf = <||>},
  Needs["FunKit`"];
  (* A large flow assembles a kernel with one summand per diagram GROUP (ZA4: 1274). Several codegen
     steps (the integrand Sum, COEN's expression lowering) recurse ~linearly in that count, so the
     default $RecursionLimit of 1024 is exceeded — and $RecursionLimit::reclim does NOT abort, it
     returns a held expression and the script continues to "DONE" having SILENTLY skipped the kernel
     (the ZA4 silent-skip). Raise the limit generously for the emission; a real runaway would still
     hit the (much higher) ceiling and surface. *)
  $RecursionLimit = Max[$RecursionLimit, 1048576];
  name  = OptionValue["Name"];  dress = OptionValue["Dressings"];
  scalarParams = OptionValue["ScalarParams"];  (* loop-independent scalar doubles threaded into the signature *)
  angleDefs = OptionValue["AngleDefs"];  crossCSE = OptionValue["CrossTraceCSE"];
  gc = OptionValue["GlobalCollect"];  (* Route-B: fold the WHOLE integrand (all diagrams, dressings as
    inert symbols) into ONE collected polynomial -> one kernel, like FORM. crossCSE is subsumed by it. *)
  decor = OptionValue["Decorator"];
  kns = OptionValue["KernelNamespace"];  sns = OptionValue["SupportNamespace"];
  runInc = OptionValue["RuntimeInclude"];  extraInc = OptionValue["ExtraIncludes"];
  interpTy = ntDressType[OptionValue["DressingType"]];
  ns    = OptionValue["Namespace"] /. Automatic -> ToLowerCase[name];
  nsHome = kns <> "::" <> ns;  (* where the generated trace fns / nenv / fill live *)
  args  = k["Args"];  frame = k["Frame"];  env = k["Env"];
  (* Vestigial compile-chain salt — see mkDirectKernel: SU(N) heads carry their own rank N. *)
  nc = 0;
  mask  = Association @ KeyValueMap[#1 -> frameMask[resolveComponents[#1, frame]] &, env];

  fillArgs   = Select[args, # =!= Global`k &];                 (* scalars the fill needs *)
  (* NUMERIC backend: build the component table over a compact parametrisation. Automatic uses the
     unit-loop spec (loop = magnitude × unit-direction symbols + unit constraint) so the contraction is
     as compact as the sp basis for BOTH Lorentz and Dirac nets. User "Components" are taken verbatim
     (polynomialised for any Sqrt/trig). *)
  ncomp = Module[{uc = OptionValue["Components"], ud = OptionValue["SymbolDefs"], pf, ad, ug},
    {pf, ad, ug} = Which[
       uc =!= Automatic, Append[polyFrameSpec[uc], {}],                       (* explicit user Components *)
       unitLoopOkQ[frame, Global`p, Global`l1], unitLoopFrameSpec[frame, Global`p, Global`l1], (* compact vacuum SP frame *)
       True, Append[polyFrameSpec[frame], {}]];                               (* finite-T / general frame *)
    symDefs = Join[ad, ud];
    numericComponents[env, pf, symDefs, ug]];
  fillArgSig = StringRiffle[("double " <> SymbolName[#]) & /@ fillArgs, ", "];

  (* walk diagrams: each one is a Lorentz trace x a colour factor x a dressing/kinematic coeff. The
     bridge distribution split each 4-gluon vertex into colour channels, so MANY diagrams share the
     SAME coeff and differ only in (colour x Lorentz). Collect per diagram {coeff, lorNet, colNet},
     then GROUP by coeff: per group the generator folds the (constant) colour into the Lorentz poly
     and combines them — combinedTr_g = Σ_d colN_d · trN_d — so the kernel evaluates ~one polynomial
     per Feynman graph (≈5), not one per channel (51). Vanishing-colour diagrams drop out for free. *)
  (* Fundamental colour (a quark loop's T^a) can't go through the adjoint SUNNet, but it IS a
     small et-instantiable contraction (unlike the four-gluon type), so we evaluate it with the et
     engine as a kernel `constexpr double` — exactly the direct path's constant-colour handling —
     and multiply it into the integrand. The generator then carries an identity colour (SUNNet{}).
     With fundamental colour present we do NOT group diagrams (each keeps its own colour constant). *)
  hasFund = ! FreeQ[k["Diagrams"], _ntSUNT | _ntSUNDeltaFund];
  (* Dirac VERTEX: projector i + imaginary non-abelian colour (f^abc T^b T^c = (iN/2)T^a). Per
     diagram coeff*colour is real (the i's combine), so we keep the colour constant COMPLEX and
     take Re of the assembled integrand — exactly the dense treatment. *)
  complexQ = ! FreeQ[k["Diagrams"], Complex];
  (* CrossTraceCSE fills a `double tarr[]` via trace_all; a complex trace would be truncated to its
     real part there, so the RePart re/im split's ntIm(tarr[i]) would read 0 and silently drop the
     imaginary contribution. GlobalCollect subsumes crossCSE, so disable it when the flow is complex. *)
  If[crossCSE && complexQ, Message[mkGenerateKernel::crosscseComplex]; crossCSE = False];
  invNets = {}; invRest = {}; colourNets = {}; colDecls = {}; colToks = {}; preamble = {};
  factorNets = {}; lorFacOf = {}; factorCompOf = <||>;  (* factor net indices, per-net factor-id list, net->factor-id *)
  $ctCache = <||>;  (* clear the compileTInv memo cache for this generation *)
  resetDiagDr[];     (* clear the per-component diagonal-dressing registry for this generation *)
  resetDr[];         (* clear the scalar-dressing (ntDressedNum) registry for this generation *)
  (* frame resolver for dressed-numerator option coefficients (compileDirac → dressedSlotStr): the same
     ntSP/ntSPS/ntVec[q,i] → component substitution used for diag["Coeff"] below. *)
  $ntDressResolve = Function[s, s /. {
     ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame],
     ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]],
     ntVec[q_, ii_Integer] :> resolveComponents[q, frame][[ii + 1]]}];
  Module[{diagData = {}},
    ntLog["[prof] per-diagram net-build (", Length[k["Diagrams"]], " diagrams): ", First@AbsoluteTiming[
    MapIndexed[Function[{diag, di}, Module[{coeff, col = "SUNNet{}", colTok = "", entries = {}, d = di[[1]] - 1,
        pureLorAcc = {}, nNonConst = 0, nNCDirCol = 0, diracComps = {},
        nFundConst = Count[diag["Components"],
          c_ /; c["Constant"] && ! FreeQ[c["Factors"], _ntSUNT | _ntSUNDeltaFund]]},
      coeff = diag["Coeff"] /. {
    ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame],
    ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]],
    ntVec[q_, i_Integer] :> resolveComponents[q, frame][[i + 1]]};
      (* one fundamental constant component -> the legacy `_col<d>` name (kernels stay byte-identical);
         several (e.g. a Yukawa loop's flavour AND colour traces) -> tag each `_col<d>_<ci>`. *)
      MapIndexed[Function[{comp, ci}, Module[{str, scal,
          tag = If[nFundConst > 1, ToString[d] <> "_" <> ToString[ci[[1]] - 1], ToString[d]]},
        If[comp["Constant"],
          (* Constant SU(N) component — colour and/or flavour, fundamental and/or adjoint: every group
             head carries its own rank N, so all fold through the SINGLE numeric SUNNet path
             (sun_value_cx, contracting each rank separately and multiplying). The fold is COMPLEX
             (-> Cx): when an imaginary non-abelian-vertex colour (f^abc T^b T^c = (iN/2) T^a) is folded
             in, the diagram's trace lowers to a complex `tr_i`; the kernel multiplies by the complex
             dressing coefficient and the consumer takes Re, so the imaginary part survives. ACCUMULATE
             components — a Yukawa loop carries SEVERAL constant components (a colour trace AND a flavour
             trace); mergeColNet folds them together so none is dropped (a bare `col = str` would keep
             only the last, e.g. the δ^ii = Nf factor, giving a ~50% wrong trace). Single-component
             flows are unchanged: mergeColNet["SUNNet{}", str] == str. *)
          ({str, scal} = compileColG[Times @@ comp["Factors"], diag["Ids"]]; coeff *= scal;
           col = mergeColNet[col, str]),
          (* Non-constant component. The DISCONNECTED components of ONE diagram MULTIPLY (their scalar
             trace values) — exactly like the dense path's `Times @@ toks` (compileDiagramDense) — they
             are NOT separate summed diagrams. Collect them so the post-loop assembly can form the
             product (see there). Route by structure:
               - any colour (entangled in a Plus, or a top-level T^a × …) or a gamma chain: collect the
                 component's splitColourGroupsInv entries (dirac_value net per colour branch / chunked
                 Lorentz net) as ONE Dirac/colour component in diracComps;
               - pure Lorentz (no colour, no gamma): accumulate factors — ALL pure-Lorentz components
                 fold into ONE product net (disjoint ids make the C++ contract_factors multiply them). *)
          (nNonConst++;
           If[colourEntangledQ[comp["Factors"]] || ! FreeQ[comp["Factors"], _ntGamma | _ntGamma5 | _ntDeltaDirac | _ntDressedNum],
             (nNCDirCol++; AppendTo[diracComps, splitColourGroupsInv[comp["Factors"], diag["Ids"], env, mask, nc]]),
             pureLorAcc = Join[pureLorAcc, comp["Factors"]]])]]],
        diag["Components"]];
      (* ---- assemble the diagram's nets from its non-constant components -----------------------------
         A diagram with K disconnected non-constant components is a PRODUCT of K independent closed
         scalars (each a Dirac/colour trace or a pure-Lorentz scalar) — exactly the dense path's
         `coeff * Times @@ toks`. ALL pure-Lorentz components fold into ONE product factor; each
         Dirac/colour component is its own factor. Two regimes:
           - <= 1 non-constant factor: the EXISTING additive path. The single Dirac component's entries
             (colour folded, GlobalCollect-fusible) OR the single combined pure-Lorentz product net is
             appended with diagData = coeff*scal. Flows without a disconnected diagram regenerate
             BYTE-IDENTICAL (any K>=2 diagram previously aborted, so none is committed).
           - >= 2 factors: FACTORED product. One Dirac component is the additive BASE (traceRef[gi],
             diagData = coeff*scal); every other component is emitted as its own fused trace GROUP
             (its per-entry scalar folded into the net, its colour folded by the group sum) and tagged
             with a factor id. The base carries the list of factor ids; the assembly multiplies in
             Π traceRef[factor groups] — P computed ONCE per component, no trace-polynomial blow-up.
         scalarleakCheck: each entry's restNet scalars (e[[4]] {restStr,scal}) become the generator's
         `dsc[]` numeric constants — a symbolic Lorentz tensor the net builder failed to fold would be
         CForm'd into undeclared C++. Catch it loudly, with the offender. *)
      Module[{nDir = Length[diracComps], hasLor = pureLorAcc =!= {}, factorComps, factorIds = {},
              scalarleakCheck, appendRec},
        scalarleakCheck = Function[es, Do[Module[{badS = FirstCase[ee[[4]], {_, s_} /; ! NumericQ[s] :> s, Missing[]]},
           If[! MissingQ[badS], Message[mkGenerateKernel::scalarleak, d, badS]; Abort[]]], {ee, es}]];
        (* rec = {colNet, cores, dData, restList, lorFac(None|{ids..}), factorId(None|id)}; returns net idx *)
        appendRec = Function[rec, Module[{ni = Length[invNets]},
           invNets = Append[invNets, rec[[2]]]; invRest = Append[invRest, rec[[4]]];
           colourNets = Append[colourNets, rec[[1]]]; colToks = Append[colToks, colTok];
           diagData = Append[diagData, rec[[3]]]; lorFacOf = Append[lorFacOf, rec[[5]]];
           If[rec[[6]] =!= None, factorNets = Append[factorNets, ni]; factorCompOf[ni] = rec[[6]]]; ni]];
        If[nDir + Boole[hasLor] <= 1,
          (* ---- single non-constant factor (or none): EXISTING additive path, byte-identical ---- *)
          Module[{baseEntries = Which[
              nDir == 1, diracComps[[1]],
              hasLor,    ({"SUNNet{}", {#[[1]]}, #[[2]], {{"", 1}}} & /@
                            chunkLorInv[Times @@ pureLorAcc, diag["Ids"], env, mask, nc]),
              True,      {}]},
            scalarleakCheck[baseEntries];
            Do[appendRec[{mergeColNet[col, e[[1]]], e[[2]], coeff e[[3]], e[[4]], None, None}],
               {e, If[baseEntries === {}, {{"SUNNet{}", {"konst(1.0)"}, 1, {{"", 1}}}}, baseEntries]}]],
          (* ---- >= 2 factors: factored product of disconnected components ----
             EVERY non-constant component becomes its OWN fused trace group: its colour-branch entries
             are summed WITHIN the group (GlobalCollect-style, colour folded by the group sum, the entry
             scalar folded into the net) so the group trace IS the component scalar. The diagram then
             contributes ONE additive ANCHOR term = coeff * colv(col) * Π(component group traces).
             Crucially this does NOT force the components' many entries into singletons (which would
             defeat colour-channel fusion and explode the trace count for a four-quark-dressed loop) —
             each component is exactly ONE trace, mirroring the dense path's per-component scalar. *)
          (factorComps = If[hasLor,
             Append[diracComps,
               {{"SUNNet{}", {#[[1]]}, #[[2]], {{"", 1}}}} &[
                  compileTInv[Times @@ pureLorAcc, diag["Ids"], env, mask, nc]]],
             diracComps];
           Do[Module[{compEntries = factorComps[[ci]], fid = Length[invNets]},
               scalarleakCheck[compEntries]; AppendTo[factorIds, fid];
               Do[appendRec[{e[[1]], e[[2]], 1, ({#[[1]], #[[2]] e[[3]]} &) /@ e[[4]], None, fid}],
                  {e, compEntries}]],
              {ci, 1, Length[factorComps]}];
           (* the anchor: a trivial unit net carrying the diagram coeff, the constant colour `col`
              (folded once, via the group sum colv(col)*1), and the list of component factor ids. *)
           appendRec[{col, {"konst(1.0)"}, coeff, {{"", 1}}, factorIds, None}])]]]], k["Diagrams"]]], " s"];
    (* ---- per-component diagonal dressings (ntSUNDiag{Fund,Adj}) ----------------------------------
       A diagram whose colour net carries a diag factor folds (via the validated C++ engine,
       sun_value_dressed, run through the build-time seam) to a SUNPoly Σ_t coeff_t Π D^{dr}_{comp}.
       The per-component dressing values D are RUNTIME (a kernel std::array indexed by the group
       component), so the colour can't be a compile-time constant baked into the trace: instead we
       (a) replace the diagram's colour net by the IDENTITY (so the generator folds colv=1 and the
       trace stays colour-free), and (b) build a runtime token `Σ_t coeff_t Π ntDiagDress[name,comp,
       scale]` multiplied into the integrand. The Dirac/Lorentz trace is still computed ONCE — the
       per-component sum is a handful of array lookups, not a diagram per component. *)
    diagTokExpr = Table[1, {Length[colourNets]}];
    dressedIdx = Select[Range[Length[colourNets]], StringContainsQ[colourNets[[#]], "SUN::diag"] &];
    If[dressedIdx =!= {},
      Module[{polys, resolveScale,
              incDir = OptionValue["IncludeDir"] /. Automatic :> resolveIncludeDir[]},
        resolveScale[s_] := s /. {
          ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame],
          ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]],
          ntVec[q_, ii_Integer] :> resolveComponents[q, frame][[ii + 1]]};
        polys = diagColPolys[colourNets[[dressedIdx]], incDir];
        MapThread[Function[{d, p},
          diagTokExpr[[d]] = Total[Function[term,
             (term[[1]] + I term[[2]]) * (Times @@ (Function[dc,
                Global`ntDiagDress[$diagDrTable[dc[[1]]]["Name"], dc[[2]],
                   resolveScale[$diagDrTable[dc[[1]]]["Scale"]]]] /@ term[[3]]))] /@ p];
          colourNets[[d]] = "SUNNet{}"], {dressedIdx, polys}]];
      ntLog["[prof] diagonal-dressed diagrams: ", Length[dressedIdx],
            " (per-component colour-sum folded via sun_value_dressed seam)"]];
    (* trace reference: with CrossTraceCSE the kernel fills a `tarr[]` once via trace_all() and reads
       tarr[i]; otherwise it calls the independent tr_i(fenv). *)
    traceRef = If[crossCSE, "tarr[" <> ToString[#] <> "]", nsHome <> "::tr" <> ToString[#] <> "(fenv)"] &;
    (* Group diagrams (0-based) into traces by colour kind. ADJOINT-colour diagrams (colTok=="" —
       their real colour was folded numerically into the generator polynomial) FUSE by identical
       dressing coeff, so the kernel evaluates ~one polynomial per Feynman graph; cross-diagram like
       terms then collect in the invariant basis (many channels collapse to a few traces this way).
       FUNDAMENTAL-colour diagrams (a complex generator T^a, kept in the kernel as _colN so the trace
       polynomial stays REAL) stay singleton, each times its own _colN. A flow that mixes both colour
       kinds thus fuses the adjoint part and keeps the fundamental part per-diagram. *)
    (* diag-dressed diagrams (diagTokExpr =!= 1) carry a per-diagram RUNTIME colour-sum token, so —
       like the fundamental-colTok diagrams — they stay singletons (they cannot fuse by dressing
       coefficient alone, the colour token differs). lorFac diagrams (lorFacOf =!= None — a Dirac trace
       times a disconnected pure-Lorentz scalar) likewise stay singletons: each carries a per-diagram
       multiplicative trace, so it must not fuse with another diagram's entries.
       The FACTOR nets (P, indices in factorNets) are EXCLUDED from the additive groups and appended as
       their own singleton trace groups at the tail of g — generated as traces but referenced only
       multiplicatively via lorFac, never summed into the integrand. nAdd marks the additive/factor
       boundary; pGroupOf maps a P net (0-based) to its (0-based) trace-group ordinal. *)
    Module[{adj, fund, additivePos, factorPos = (# + 1) & /@ factorNets, gAdd, gFactor},
      additivePos = Complement[Range[Length[diagData]], factorPos];
      adj  = Select[additivePos, colToks[[#]] === "" && diagTokExpr[[#]] === 1 && lorFacOf[[#]] === None &];
      fund = Select[additivePos, colToks[[#]] =!= "" || diagTokExpr[[#]] =!= 1 || lorFacOf[[#]] =!= None &];
      gAdd = Join[(# - 1) & /@ GatherBy[adj, diagData[[#]] &], List /@ (fund - 1)];
      (* each disconnected factor COMPONENT (possibly several colour-branch nets) fuses into ONE trace
         group, so its group trace = the component scalar (colour folded by the group sum). *)
      gFactor = GatherBy[factorNets, factorCompOf[#] &];   (* factorNets are 0-based net indices *)
      g = Join[gAdd, gFactor]; nAdd = Length[gAdd];
      pGroupOf = Association[MapIndexed[
         (factorCompOf[#1[[1]]] -> (nAdd + #2[[1]] - 1)) &, gFactor]]];
    (* GlobalCollect folds ALL colour numerically (colToks all empty above), so this groups EVERY
       diagram by its dressing coefficient — the quark-loop colour channels that the legacy path kept
       as fundamental singletons now MERGE into their Feynman graph, collapsing the trace count toward
       FORM's handful. The dressing coefficient stays FACTORED in `diagData` (COEN CSEs it, like FORM's
       _repl), so each group is one collected kinematic trace × its dressing — not a flat polynomial. *)
    (* Sum the ADDITIVE groups only (1..nAdd); factor groups (nAdd+1..) are referenced multiplicatively
       via lorFac. A factored (disconnected) diagram multiplies in its other components' scalars
       Π traceRef[factor groups] (each computed ONCE, as a separate trace), exactly like the dense
       `coeff * Times @@ component-scalars`. *)
    integrand = Sum[With[{rep = g[[gi, 1]]},
        diagData[[rep + 1]] * If[colToks[[rep + 1]] === "", 1, colToks[[rep + 1]]] *
          diagTokExpr[[rep + 1]] *
          If[lorFacOf[[rep + 1]] === None, 1,
             Times @@ (traceRef[pGroupOf[#]] & /@ lorFacOf[[rep + 1]])] *
          traceRef[gi - 1]],
      {gi, nAdd}]];
    (* OPTIMIZATION (opt-in, "PruneRealTraces"): a group whose dressing coefficient is REAL has only
       Re(trace) consumed (the consumer takes Re of the whole kernel; a real coeff cannot move
       Im(trace) into the real part). Flag it so the generator emits a `double` trace and never
       computes the dead imaginary half. Empty list (default) => generator emits all-complex traces,
       which the probe needs to verify cancellation — see the "PruneRealTraces" note. *)
    realOnlyG = If[TrueQ[OptionValue["PruneRealTraces"]],
       (FreeQ[diagData[[#[[1]] + 1]], Complex]) & /@ g, {}];
  nGrp = Length[g];
  (* [B2 localize] surface the post-net-build shape so a silent no-output (e.g. empty nets / empty
     grouping) is visible rather than appearing as a clean DONE. *)
  ntLog["[prof] post-net-build: nets=", Length[invNets], " groups(nGrp)=", nGrp, " complexQ=", complexQ];
  If[Length[invNets] === 0 || nGrp === 0,
    Message[mkGenerateKernel::emptynets, name, Length[invNets], nGrp]; Abort[]];

  (* kinematic angle defs (kept symbolic in the dressing): emit once as named temporaries. *)
  angleDecls = ("const double " <> SymbolName[First[#]] <> " = " <> cppFlat[Last[#]] <> ";") & /@ angleDefs;
  (* NB: deliberately NO `using std::complex;` — unqualified complex<double> resolves to the
     support namespace's `complex` alias. That indirection lets a device/CUDA support header
     substitute a device-safe complex (std::complex arithmetic lowers to gcc _Complex builtins
     that nvcc silently miscompiles to 0 in device code), without changing the emitted kernel. *)
  (* the fenv setup block: declare fenv, (dressed only) compute each dressing atom into dr_<id>, fill,
     and (CrossTraceCSE) precompute the traces. *)
  With[{hasDr = ! FreeQ[invNets, _ntDressedCore]},
   Module[{coreBlock},
    coreBlock = {"double fenv[(" <> nsHome <> "::nenv) > 0 ? (" <> nsHome <> "::nenv) : 1];",
       Sequence @@ If[hasDr,
         KeyValueMap[Function[{id, atom},
           "const double dr_" <> ToString[id] <> " = " <> cppFlat[atom] <> ";"], $drTable], {}],
       With[{fillCallArgs = If[hasDr,
              Join[SymbolName /@ fillArgs, ("dr_" <> ToString[#]) & /@ Sort[Keys[$drTable]]],
              SymbolName /@ fillArgs]},
         nsHome <> "::fill(fenv, " <> StringRiffle[fillCallArgs, ", "] <> ");"],
       If[crossCSE,
         "double tarr[" <> ToString[nGrp] <> "]; " <> nsHome <> "::trace_all(fenv, tarr);", Nothing]};
    (* DRESSED: the dr_<id> dressing expressions can reference the derived kinematic angles, so the
       angle (and colour) decls must precede the fenv block. NON-dressed: keep the original order
       (fenv before angle/colour decls) so those kernels regenerate byte-identical. *)
    preamble = StringRiffle[If[hasDr,
       Join[ntSupportUsings[sns], angleDecls, colDecls, coreBlock, preamble],
       Join[ntSupportUsings[sns], coreBlock, angleDecls, colDecls, preamble]], "\n"]]];

  mkParam[nm_, ty_] := <|"Name" -> If[StringQ[nm], nm, SymbolName[nm]], "Type" -> ty, "Const" -> True, "Reference" -> True|>;
  (* a PER-COMPONENT dressing (ntSUNDiag{Fund,Adj}) is passed as a std::array<interp, dim> indexed by
     the group component; a plain dressing keeps the scalar interpolator type. *)
  diagDressDims = Association[(#["Name"] -> #["Dim"]) & /@ Values[$diagDrTable]];
  With[{dressTy = Function[nm, If[KeyExistsQ[diagDressDims, nm],
      "std::array<" <> interpTy <> ", " <> ToString[diagDressDims[nm]] <> ">", interpTy]]},
    kernelParams = Join[mkParam[#, "double"] & /@ args, mkParam[#, "double"] & /@ scalarParams, mkParam[#, dressTy[#]] & /@ dress];
    constParams  = Join[mkParam[#, "double"] & /@ Select[args, # === Global`p || # === Global`k &],
                        mkParam[#, "double"] & /@ scalarParams,
                        mkParam[#, dressTy[#]] & /@ dress];
    (* dressed kernels: fill() takes one `double dr_<id>` per dressing atom — the kernel body computes
       the atom's value (regulators / interpolators in scope there) and passes it. Matches fm.dress. *)
    If[! FreeQ[invNets, _ntDressedCore],
       fillArgSig = fillArgSig <> StringJoin[(", double dr_" <> ToString[#]) & /@ Sort[Keys[$drTable]]]]];

  kernelFn = FunKit`MakeCppFunction[integrand, "Name" -> "kernel", "Prefix" -> decor,
    "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> kernelParams, "Body" -> preamble];
  constFn = ntConstFn[OptionValue["Constant"], decor, constParams, sns];
  classStr = FunKit`MakeCppClass["TemplateTypes" -> {"REG"}, "Name" -> name,
    "MembersPublic" -> If[TrueQ[OptionValue["RegulatorAlias"]],
      Prepend[{kernelFn, constFn}, "using Regulator = REG;"], {kernelFn, constFn}],
    "MembersPrivate" -> {privDefs[decor]}];
  hdrInc = FileNameTake[headerFile];
  header = FunKit`MakeCppHeader[
    (* the numeric kernel is flat straight-line arithmetic: the generated trace functions (hdrInc) plus
       the support runtime; no tensor-engine headers. *)
    "Includes" -> Join[extraInc, ntRuntimeIncludes[runInc],
       {"numtracer/sun/sun_data.hpp", "numtracer/expr/real_cse.hpp", hdrInc}],
    "Body" -> ntWrapBody[kns, classStr, name]];

  (* emit the generator source (the numeric matrix-product backend is the single generation path). *)
  ntLog["[prof] emitNumericGenerator: ", First@AbsoluteTiming[
  {genPre, genUnits, genDecl, genMain} =
    emitNumericGenerator[invNets, invRest, colourNets, g, ncomp, ns, fillArgSig, kns, complexQ, realOnlyG];], " s"];
  (* Split generator: a main TU + N net-builder unit TUs + a decl header (all in the tests/gen/ dir), so the
     net-builder codegen compiles in parallel (see emitNumericGenerator). The main `#include`s the decl. *)
  declFile = StringReplace[genFile, ".cpp" -> "_nets.hh"];
  unitFiles = Table[StringReplace[genFile, ".cpp" -> "_u" <> ToString[u - 1] <> ".cpp"], {u, 1, Length[genUnits]}];
  Export[declFile, genDecl, "Text"];
  (* each unit #includes the shared decl header so its net builders can call the cross-unit CSE
     accessors (lc<k>()/dc<k>()) and sibling net builders, with the declarations parsed once per TU. *)
  Module[{uInc = "#include \"" <> FileNameTake[declFile] <> "\"\n"},
    Do[Export[unitFiles[[u]], uInc <> genUnits[[u]], "Text"], {u, 1, Length[genUnits]}]];
  genSrc = genPre <> "\n#include \"" <> FileNameTake[declFile] <> "\"\n\n" <> genMain;
  Export[genFile, genSrc, "Text"];
  Print["wrote generator: ", genFile, " (+ ", Length[genUnits], " net units + decl header)"];

  (* run the generator at codegen time -> the committed straight-line kernel header. The binary's
     stdout is redirected straight to the header FILE via the shell (Run), not captured in memory
     by RunProcess — the A4 kernel header is ~40k lines, and in-memory capture is fragile at that
     size. *)
  If[OptionValue["RunGenerator"],
    incDir = OptionValue["IncludeDir"] /. Automatic :> resolveIncludeDir[];
    bin = FileNameJoin[{$TemporaryDirectory, "gen_" <> ns}];
    (* Time COMPILE and RUN separately — both count toward generation time, but the levers differ
       (compile: TU size / templates; run: reduce+rebase). Reported so neither is hidden. Compile the
       main TU (-O2) and the net-builder units (-O0) CONCURRENTLY (`&` + `wait`), then link — the unit
       codegen runs across cores. A failed unit compile leaves its .o missing -> the link rc is nonzero,
       so the existing rc check catches it. *)
    Module[{tcc, cc, mainObj, unitObjs, pcmd, lcmd, clog = bin <> "_compile.log", cxx = resolveGenCxx[],
            mainOpt},
      mainObj = bin <> "_main.o";
      unitObjs = Table[bin <> "_u" <> ToString[u - 1] <> ".o", {u, 1, Length[unitFiles]}];
      (* Main-TU optimisation level. The generator runs ONCE (~1 s), so -O2 on it is wasteful — but its
         contraction loop must stay inlined or the RUN slows ~10x. For the DRESSED (collected) path the
         main TU is heavy to compile at -O2 (the DPoly templates + the inline dch/dsl chain/slot literals
         instantiate slowly: measured 28.7 s at -O2 vs 17.5 s at -O1 with an IDENTICAL ~1 s run, vs -O0
         which compiles in 6 s but slows the run to ~9 s). -O1 keeps the contraction optimised while
         dropping -O2's costly extra passes, so it is the dressed default. Non-dressed stays -O2
         (already fast to compile, and its reduce+rebase RUN on big flows benefits from -O2). Override
         with NT_GEN_MAIN_OPT. *)
      mainOpt = With[{e = Environment["NT_GEN_MAIN_OPT"]},
        Which[StringQ[e] && e =!= "", e,
          ! FreeQ[invNets, _ntDressedCore], "-O1",
          True, "-O2"]];
      (* RAM-bounded parallel compile: run at most $ntCompileJobs `cxx` at once (xargs -P), and cap
         each at ~17 GB virtual (ulimit -v). Peak RAM ~ jobs x per-unit; with the unit count scaled so
         each TU is small (~12 net-builders) this stays well under the machine limit for any flow size.
         Both phases redirect compiler stdout+stderr to `clog` (compile truncates, link appends) so the
         genfail message can quote the actual g++ diagnostic, not just the exit code. *)
      pcmd = "printf '%s\\0' " <> StringRiffle[("\"" <> # <> "\"") & /@ Join[
          {"(ulimit -v 17000000; " <> cxx <> " -std=c++20 -ftemplate-depth=4000 " <> mainOpt <> " -pthread -I '" <> incDir <> "' -c '" <> genFile <> "' -o '" <> mainObj <> "')"},
          Table["(ulimit -v 17000000; " <> cxx <> " -std=c++20 -ftemplate-depth=4000 -O0 -I '" <> incDir <> "' -c '" <> unitFiles[[u]] <> "' -o '" <> unitObjs[[u]] <> "')",
            {u, 1, Length[unitFiles]}]], " "] <> " | xargs -0 -P " <> ToString[$ntCompileJobs] <> " -I CMD bash -c CMD > '" <> clog <> "' 2>&1";
      lcmd = cxx <> " -pthread '" <> mainObj <> "' " <> StringRiffle[("'" <> # <> "'") & /@ unitObjs, " "] <> " -o '" <> bin <> "' >> '" <> clog <> "' 2>&1";
      {tcc, cc} = AbsoluteTiming[Run[pcmd]; Run[lcmd]];
      ntLog["[time]   generator compile (", cxx, ", ", Length[unitFiles], " parallel units + main): ", tcc, " s"];
      If[cc =!= 0,
        Module[{lines = If[FileExistsQ[clog], StringSplit[ReadString[clog], "\n"], {}], head, nshow = 50},
          head = Take[lines, UpTo[nshow]];
          Message[mkGenerateKernel::genfail, cxx <> " compile/link rc=" <> ToString[cc] <>
            "\n--- first " <> ToString[Length[head]] <> " of " <> ToString[Length[lines]] <> " line(s) of " <> clog <> " ---\n" <>
            StringRiffle[head, "\n"] <>
            If[Length[lines] > nshow, "\n... (" <> ToString[Length[lines] - nshow] <> " more line(s) in " <> clog <> ")", ""]];
          Abort[]]]];
    (* run into a TEMP file, validate (rc==0 AND non-empty), then move into place — so a crashed
       generator (e.g. thread-limited Run[]) never silently truncates the committed header. *)
    Module[{tmp = headerFile <> ".tmp", rc, sz, trun},
      {trun, rc} = AbsoluteTiming[
        Run["'" <> bin <> "' -n '" <> ns <> "' -d '" <> decor <> "'" <>
          If[OptionValue["FullParallel"], " -p", ""] <> " > '" <> tmp <> "'"]];
      ntLog["[time]   generator run (reduce+rebase+lower): ", trun, " s"];
      sz = If[FileExistsQ[tmp], FileByteCount[tmp], 0];
      If[rc =!= 0 || sz < 64,
        If[FileExistsQ[tmp], DeleteFile[tmp]];
        Message[mkGenerateKernel::genfail, "generator run rc=" <> ToString[rc] <> " bytes=" <> ToString[sz] <>
          " (committed header left intact)"]; Abort[]];
      CopyFile[tmp, headerFile, OverwriteTarget -> True]; DeleteFile[tmp];
      Print["wrote header: ", headerFile, " (", sz, " bytes)"]]];

  (* semantic complexQ: with the real traces now generated, probe whether Im(integrand) actually
     vanishes (the projector-i / colour-i factors usually cancel to a real flow). If so, re-emit a
     REAL (double) kernel — losslessly, since Im≡0 — so the real DiFfRG integrators bind directly and
     no std::complex survives into device code. Only meaningful once the traces exist (RunGenerator). *)
  If[complexQ && TrueQ[OptionValue["RunGenerator"]] && TrueQ[OptionValue["RealProbe"]],
    Module[{verdict = numericImagProbeRealQ[integrand, args, fillArgs, angleDefs, angleDecls, nsHome, headerFile, $drTable],
            extraPriv = {}},
      Switch[verdict,
        "Pure",   (* real value: project to Re. Wrap trace tokens in ntRe FIRST so the kernel is
                     provably double-typed even if a trace function is complex-typed — the probe
                     proved |full-proj|≈0, i.e. Σ Im(c)·tr ≈ 0, so dropping the ntIm terms is exact
                     and ntRe(double)=passthrough / ntRe(complex)=.real() keeps it real arithmetic. *)
          ntLog["[probe] real & projection exact -> clean REAL (double) kernel"];
          integrand = (integrand /. s_String :> Global`ntRe[s]) /. Complex[a_, b_] :> a;
          complexQ = False;  extraPriv = {ntReImDefs[decor]},
        "RePart", (* real value but a complex trace: re/im split -> Σ[Re(c)·tr.real() − Im(c)·tr.imag()] *)
          ntLog["[probe] real value via complex trace(s) -> double kernel (re/im split)"];
          (* The integrand is LINEAR in the trace tokens (strings). Re(Σ c·tr) = Σ[Re(c)·ntRe(tr) −
             Im(c)·ntIm(tr)]. The naive `(… /. s:>ntRe[s]+I ntIm[s]) /. Complex[a_,b_]:>a` is WRONG:
             Mathematica keeps `i·X·(ntRe+I ntIm)` as an UNEXPANDED product, so `/. Complex:>a` zeroes
             the leading `i` factor and DROPS the whole term — silently losing the real −X·ntIm(tr)
             contribution of every complex trace. Instead split each token's coefficient into real /
             imaginary parts via an `ii`-substitution (I → real symbol ii), which keeps the dressing
             coefficients FACTORED (unlike ComplexExpand, which un-factors them and defeats COEN's CSE). *)
          integrand = Module[{toks = Union[Cases[integrand, _String, Infinity]], tsym, lin, ii},
            tsym = AssociationThread[toks -> Table[Unique["tr$"], {Length[toks]}]];
            lin = integrand /. tsym;
            Total[Function[t, Module[{cc = Coefficient[lin, tsym[t]] /. Complex[ar_, ai_] :> ar + ii*ai},
               Coefficient[cc, ii, 0]*Global`ntRe[t] - Coefficient[cc, ii, 1]*Global`ntIm[t]]] /@ toks]
              + ((lin /. Thread[Values[tsym] -> 0]) /. Complex[ar_, ai_] :> ar)];
          complexQ = False;  extraPriv = {ntReImDefs[decor]},
        _,
          ntLog["[probe] imaginary part survives -> keeping the complex kernel (consumer takes Re)"]];
      If[verdict === "Pure" || verdict === "RePart",
        kernelFn = FunKit`MakeCppFunction[integrand, "Name" -> "kernel", "Prefix" -> decor,
          "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> kernelParams, "Body" -> preamble];
        classStr = FunKit`MakeCppClass["TemplateTypes" -> {"REG"}, "Name" -> name,
          "MembersPublic" -> If[TrueQ[OptionValue["RegulatorAlias"]],
            Prepend[{kernelFn, constFn}, "using Regulator = REG;"], {kernelFn, constFn}],
          "MembersPrivate" -> Join[{privDefs[decor]}, extraPriv]];
        header = FunKit`MakeCppHeader[
          "Includes" -> Join[extraInc, ntRuntimeIncludes[runInc],
             {"numtracer/sun/sun_data.hpp", "numtracer/expr/real_cse.hpp", hdrInc}],
          "Body" -> ntWrapBody[kns, classStr, name]]]]];

  (* kernel header (write-if-changed). *)
  If[FileExistsQ[kernelFile] && Import[kernelFile, "Text"] === header,
    Print["unchanged: ", kernelFile],
    Export[kernelFile, header, "Text"]; Print["wrote kernel: ", kernelFile]];
  kernelFile
];

(* ---- MakeNTKernel: the single public kernel emitter. ---------------------------------
   "Backend" selects the path:
     "Numeric"  — the matrix-product backend, the single canonical generation path: a build-time
                  generator program (genFile, + net-builder units + decl header), run to produce the
                  committed straight-line traces header (tracesFile), and the kernel header
                  (kernelFile) that fills the fundamental symbols and calls the traces;
     "Dense"    — brute-force numeric DTensor oracle (the entry-for-entry cross-check).
   "Dense" emits ONE file; "Numeric" emits THREE. "Backend" -> Automatic (the default) resolves
   from the argument count (one file -> Dense, three -> Numeric). The retired "ET" (symbolic) and
   "Generate" (invariant-basis) backends now error. All other options are forwarded to the selected
   path (see Options[mkDirectKernel] / Options[mkGenerateKernel] for the per-path sets). *)

Options[MakeNTKernel] = {"Name" -> "nt_kernel", "Namespace" -> Automatic,
  "Dressings" -> {}, "ScalarParams" -> {}, "Backend" -> Automatic,
  "Decorator" -> "static inline", "IncludeDir" -> Automatic, "RunGenerator" -> True,
  "FullParallel" -> False, "AngleDefs" -> {}, "CrossTraceCSE" -> False, "GlobalCollect" -> True,
  "NumericContract" -> False, "Components" -> Automatic, "SymbolDefs" -> <||>,
  "RuntimeInclude" -> "numtracer/codegen/runtime.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "numtracer_kernels", "SupportNamespace" -> "numtracer",
  "DressingType" -> Automatic, "RegulatorAlias" -> False, "RealProbe" -> True,
  "PruneRealTraces" -> False, "Constant" -> 0.};

MakeNTKernel::disconnectmix = "Diagram `1` disconnects into >= 2 Dirac/colour trace components (a \
product of independent Dirac traces, a genuine >=2-loop structure). The numeric backend handles a \
single Dirac/colour trace times any number of disconnected pure-Lorentz scalars (factored, like the \
dense backend), but does not yet multiply two or more independent Dirac traces; the dense backend does.";
MakeNTKernel::backend = "Unknown \"Backend\" `1`: use \"Numeric\" (production) or \"Dense\" (oracle).";
MakeNTKernel::nfiles = "\"Backend\" -> \"`1`\" emits `2`: MakeNTKernel[ntk, `3`, ...].";
MakeNTKernel::retired = "\"Backend\" -> \"Generate\" (the symbolic invariant-basis path) is RETIRED; \
use \"Backend\" -> \"Numeric\" (the matrix-product backend, the single canonical generation path).";
MakeNTKernel::etretired = "\"Backend\" -> \"ET\" (the symbolic expression-tensor path) has been REMOVED; \
use \"Backend\" -> \"Numeric\" (production) or \"Backend\" -> \"Dense\" (the numeric oracle).";

(* Forward the RESOLVED options (explicit opts first, then MakeNTKernel's own defaults) to the
   selected worker. Explicit opts win (first match); the defaults carry through anything not
   passed — so `SetOptions[MakeNTKernel, ...]` (used by the test setups to opt into the shim +
   the DiFfRG namespace) propagates to the workers without editing every call site. *)
MakeNTKernel[ntk : NTKernel[_], file_, opts : OptionsPattern[]] := Module[
  {be = OptionValue["Backend"] /. Automatic -> "Dense"},
  Which[
    be === "Dense",
      mkDirectKernel[ntk, file, "Backend" -> be,
        Sequence @@ FilterRules[Join[{opts}, Options[MakeNTKernel]], Options[mkDirectKernel]]],
    be === "Numeric",
      Message[MakeNTKernel::nfiles, be, "three files", "genFile, kernelFile, tracesFile"]; Abort[],
    be === "Generate",
      Message[MakeNTKernel::nfiles, be, "three files", "genFile, kernelFile, tracesFile"]; Abort[],
    be === "ET", Message[MakeNTKernel::etretired]; Abort[],
    True, Message[MakeNTKernel::backend, be]; Abort[]]];

MakeNTKernel[ntk : NTKernel[_], genFile_, kernelFile_, tracesFile_, opts : OptionsPattern[]] := Module[
  {be = OptionValue["Backend"] /. Automatic -> "Numeric"},
  Which[
    be === "Numeric",
      mkGenerateKernel[ntk, genFile, kernelFile, tracesFile, "NumericContract" -> True,
        Sequence @@ FilterRules[Join[{opts}, Options[MakeNTKernel]], Options[mkGenerateKernel]]],
    be === "Generate",
      Message[MakeNTKernel::retired]; Abort[],
    be === "ET", Message[MakeNTKernel::etretired]; Abort[],
    be === "Dense",
      Message[MakeNTKernel::nfiles, be, "one file", "file"]; Abort[],
    True, Message[MakeNTKernel::backend, be]; Abort[]]];
