# step-08: The DSL, head by head

*Builds on: [step-06](step-06.md) · Built on by: [step-11](step-11.md), [step-19](step-19.md) ·
Tags: `dsl`, `reference` · **Tier B** (a Wolfram kernel and FunKit)*

## Introduction

This step is the vocabulary reference — but a *runnable* one. Every head gets a network small
enough that its value is obvious, a generated kernel, and a check against the closed form. Reading
a table of heads teaches you their names; watching each one produce the number you predicted
teaches you what they mean.

Everything below shares the one-angle frame from [step-07](step-07.md),

$$
q_1 = (p,0,0,0), \qquad q_l = (l_1 c,\; l_1 s,\; 0,\; 0), \qquad c=\cos_1,\ s=\sqrt{1-c^2},
$$

so the only invariants available are $q_1\!\cdot\!q_1 = p^2$, $q_l\!\cdot\!q_l = l_1^2$, and
$q_1\!\cdot\!q_l = p\,l_1 c$. Every expected value is written in those three and nothing else.

### The two kinds of factor

Before the table, the distinction that governs everything: a network factor is either a **tensor**
or a **scalar**.

* A **tensor head** carries index labels and participates in contraction. `ntMetric`, `ntVec`,
  the projectors, the gammas, the SU($N$) heads.
* A **scalar** carries no labels. It multiplies a diagram and is emitted as ordinary arithmetic
  outside the trace. `ntSP`, `ntSPS`, dressings, plain numbers.

You saw the consequence in [step-06](step-06.md): the $1/p^2$ came out as `powr<-2>(p)` *multiplying*
`tr0(fenv)`, never entering it. That is the
$\sum_{\text{diagrams}}\text{coeff}\times\text{trace}$ shape again.

### The complete vocabulary

| Head | Meaning | Where it is used |
|---|---|---|
| **Lorentz** | | |
| `ntMetric[mu, nu]` | Euclidean $\delta_{\mu\nu}$ | below |
| `ntVec[q, mu]` | component $q_\mu$ (symbolic `mu`) | below |
| `ntVec[q, i]` | **integer** `i`: the *scalar* component $q_i$ | [step-16](step-16.md) |
| `ntTransProj[q, mu, nu]` | $P^T_{\mu\nu}(q) = \delta_{\mu\nu} - q_\mu q_\nu/q^2$ | below |
| `ntLongProj[q, mu, nu]` | $P^L_{\mu\nu}(q) = q_\mu q_\nu/q^2$ | below |
| `ntElectricProj[q, mu, nu]` | finite-$T$ $P^E = P^T - P^M$ | [step-16](step-16.md) |
| `ntMagneticProj[q, mu, nu]` | finite-$T$ spatial-transverse $P^M$ | [step-16](step-16.md) |
| `ntUnitVec[i]` | constant unit basis vector $e_i$ — **engine-introduced**, not written by hand | [step-19](step-19.md) |
| `ntSpatialVec[q]` | the *spatial part* of $q$ as a momentum, $\{0, q_1, q_2, q_3\}$ — **engine-introduced** from FormTracer's `vecs` | [step-16](step-16.md) |
| **Dirac** | | |
| `ntGamma[mu, din, dout]` | $\gamma^\mu$ with a Lorentz axis and two spinor axes | below |
| `ntGamma5[din, dout]` | $\gamma_5$ | extensions |
| `ntDeltaDirac[din, dout]` | spinor identity; closes a trace | below |
| `ntSigma[legA, legB, din, dout]` | the bare $\gamma$-commutator $[A,B]$ — **internal**, folded from a $\sigma^{\mu\nu}$ vertex | [step-14](step-14.md) |
| **SU($N$)** (rank is the **first** argument) | | |
| `ntSUNT[N, a, i, j]` | fundamental generator $(T^a)_{ij}$ | below |
| `ntSUNf[N, a, b, c]` | structure constant $f^{abc}$ | below |
| `ntSUNDeltaFund[N, i, j]` | fundamental $\delta_{ij}$ | extensions |
| `ntSUNDeltaAdj[N, a, b]` | adjoint $\delta^{ab}$ | extensions |
| `ntSUNDiagFund[N, i, j, spec]` | fundamental $\delta$ with a **per-component dressing** | [step-17](step-17.md) |
| `ntSUNDiagAdj[N, a, b, spec]` | adjoint $\delta$ with a per-component dressing | [step-17](step-17.md) |
| `ntEpsFund[N, i1, …, iN]` | fundamental Levi-Civita; contracted in **pairs** into $\delta$ determinants, never reaches C++ | extensions |
| **Scalars** | | |
| `ntSP[a, b]` | scalar product $a\cdot b$ | below |
| `ntSPS[a, b]` | *spatial* scalar product $\vec a\cdot\vec b$ | [step-16](step-16.md) |
| `f[args]` + `"Dressings" -> {f}` | an opaque callable, emitted as `f(args)` | below |

```{admonition} The SU(N) rank is the first argument, and it is not optional
:class: important
`ntSUNf[3, a, b, c]`, not `ntSUNf[a, b, c]`. Carrying the rank per-head is what lets a colour SU(3)
sector and a flavour SU(2) sector coexist in one network — they select different typed-out tables.
`NumTrace` validates every rank up front and aborts with a clear message if one is not a
compile-time integer, which is a much better failure than baking an unbound symbol into a kernel.
```

## The commented program

`Tutorials/step-08-dsl-heads/dsl_heads.wls`. Each `emit` call generates one kernel; the shared
helper asserts a nonzero diagram count first.

### Lorentz

```{literalinclude} ../../../Tutorials/step-08-dsl-heads/dsl_heads.wls
:language: mathematica
:start-after: "@snip begin: lorentz"
:end-before: "@snip end: lorentz"
```

### Dirac

```{literalinclude} ../../../Tutorials/step-08-dsl-heads/dsl_heads.wls
:language: mathematica
:start-after: "@snip begin: dirac"
:end-before: "@snip end: dirac"
```

Note how a spinor trace is written: the spinor labels form a **cycle**. In
`ntGamma[a, d1, d2] ntGamma[b, d2, d1]` the labels go `d1 → d2` then `d2 → d1`, closing the loop.
Get that cycle wrong — say `d2 → d3` with `d3` never used again — and you have an open spinor index,
which is not a trace at all.

### SU($N$)

```{literalinclude} ../../../Tutorials/step-08-dsl-heads/dsl_heads.wls
:language: mathematica
:start-after: "@snip begin: sun"
:end-before: "@snip end: sun"
```

Both colour nets are multiplied by `ntSP[q1, q1]`. That is not physics — it is so the kernel is not
a compile-time constant, and therefore actually exercises `fill` and the env layout.

### Scalars and dressings

```{literalinclude} ../../../Tutorials/step-08-dsl-heads/dsl_heads.wls
:language: mathematica
:start-after: "@snip begin: scalars"
:end-before: "@snip end: scalars"
```

A dressing is an **ordinary Mathematica head applied to arguments** — `myZ[p]` — that you also name
in `"Dressings"`. That declaration is the entire mechanism: it turns `myZ` into a kernel parameter
and emits the call site as `myZ(p)`.

````{admonition} `ntDress` does not currently work — use the plain-symbol form
:class: warning
`NumTracer.m` documents a head `ntDress[head, args…]` as "an opaque scalar dressing emitted
verbatim as `head(args)`". It is not implemented: no rule in `DSL.m` or `Codegen*.m` rewrites it, so
it reaches the emitter intact and produces a literal call

```cpp
const auto _interp1 = ntDress(myZ, p);   // 'ntDress' was not declared in this scope
```

which no runtime header declares. Write `myZ[p]` and declare `myZ` in `"Dressings"` instead — the
form every production flow uses.
````

The generated signature is worth seeing:

```cpp
static inline auto kernel(const double& l1, const double& cos1, const double& p, const auto& myZ)
{
  …
  const auto _interp1 = myZ(p);
  return _interp1 * tr0(fenv) * cos1 * l1 * p;
}
```

`const auto&` — the dressing is a *template* parameter, so anything callable works: a lambda, a
spline interpolator, a GPU-resident texture. The engine has no idea which, and never evaluates it
symbolically. Note also that `myZ(p)` was hoisted into `_interp1` and computed once; a dressing
called at the same argument in twenty diagrams is evaluated once per kernel call.

## Results

```bash
cmake --build build --target dsl_heads && ./build/dsl_heads
```

```text
DSL heads at (p, l1, cos1) = (1.7, 0.9, 0.35)
Lorentz
  ntMetric       =         2.89   expected         2.89   q1.q1                        ok
  ntTransProj    =      2.53598   expected      2.53598   q1.q1 - (q1.ql)^2/ql.ql      ok
  ntLongProj     =     0.354025   expected     0.354025   (q1.ql)^2 / ql.ql            ok
    P^T + P^L    =         2.89   expected         2.89   == the metric case           ok
Dirac
  ntDeltaDirac   =            4   expected            4   tr(1) = 4                    ok
  ntGamma x2     =        2.142   expected        2.142   4 (q1.ql)                    ok
  ntGamma x4     =     -7.06952   expected     -7.06952   4[2(q1.ql)^2 - q1^2 ql^2]    ok
SU(N)
  ntSUNT         =        11.56   expected        11.56   tr(T^a T^a) q1^2 = 4 p^2     ok
  ntSUNf         =        69.36   expected        69.36   f^abc f^abc q1^2 = 24 p^2    ok
scalars and dressings
  dressing       =      3.78387   expected      3.78387   myZ(p) (q1.ql) (q1.q1)       ok
ALL TESTS PASSED
```

Three of these deserve comment.

**`P^T + P^L` = the metric case.** This is the only line that does not compare against a
closed form — it compares three *independently generated kernels* against each other. Identities
between kernels are the strongest checks available, because they cannot be fooled by a sign
convention you got wrong in both the code and the oracle. [step-11](step-11.md) and
[step-22](step-22.md) build on exactly this idea.

**`tr(1) = 4`.** A constant, and the most boring line here — but a wrong $\mathrm{tr}\,\mathbb{1}$
silently rescales every fermionic flow by a power of 4, and a collapsed Dirac loop that should have
given $\mathrm{tr}\,\mathbb{1}$ and instead gave 1 is a real failure mode this engine has had. A
ratio of exactly 0.25 in a flow is its signature.

**The four-gamma value is negative.** $4[2(q_1\!\cdot\!q_l)^2 - q_1^2 q_l^2]$ at these kinematics
is $-7.07$. There is nothing positive-definite about a Dirac trace, and a sign that "looks wrong"
usually is not.

## Possibilities for extensions

1. **$\gamma_5$ and the plane.** Add
   `ntGamma5[d1,d2] ntGamma[a,d2,d3] ntGamma[b,d3,d4] ntGamma[c,d4,d5] ntGamma[e,d5,d1]` with four
   momenta. In `propFrame` every momentum lies in a plane, so the $\varepsilon$-tensor contraction
   vanishes and `NumTrace` reports **0 diagrams** — the assertion in the emit helper will fire.
   Move to `sp4Frame` and watch it become nonzero. This is the cleanest demonstration that
   "0 diagrams" is information, not a bug.

2. **The fundamental deltas.** Emit `ntSUNDeltaFund[3,i,j] ntSUNDeltaFund[3,j,i]` (expect 3) and
   `ntSUNDeltaAdj[3,a,b] ntSUNDeltaAdj[3,b,a]` (expect 8). Then check
   $T^a_{ij}\delta_{ji} = \mathrm{tr}\,T^a = 0$ — a head that gives zero, and a good test of your
   nerve when a kernel legitimately vanishes.

3. **Two groups at once.** Put an SU(3) colour factor and an SU(2) flavour factor in the same
   network — `ntSUNT[3,a,ii,jj] ntSUNT[3,a,jj,ii] ntSUNDeltaFund[2,f1,f2] ntSUNDeltaFund[2,f2,f1]`.
   Expect $4\times 2 = 8$. Then give the flavour deltas the *same labels* as the colour ones and
   reason about what happens before running it.

4. **Levi-Civita in pairs.** `ntEpsFund[3,i,j,k] ntEpsFund[3,i,j,k]` should give $3! = 6$. Now try a
   *single* $\varepsilon$ with free indices — a lone epsilon is not an SU($N$) invariant, and the
   front-end expands only pairs, so this is a case where the DSL will not let you write something
   meaningless.

5. **Break it deliberately — a broken spinor cycle.** Change `ntGamma[b, d2, d1]` to
   `ntGamma[b, d2, d3]`. The spinor line no longer closes. Read the message and compare it with the
   Lorentz-side "unbalanced label" failure from [step-06](step-06.md) extension 5.

6. **Break it deliberately — a missing rank.** Write `ntSUNf[Nc, a, b, c]` with `Nc` left as an
   unbound symbol. `NumTrace` should abort with a specific message rather than emitting a kernel.
   Confirm that it does; this guard is the reason the rank is a required leading argument.

## The plain program

```{literalinclude} ../../../Tutorials/step-08-dsl-heads/dsl_heads.wls
:language: mathematica
```

```{literalinclude} ../../../Tutorials/step-08-dsl-heads/check_dsl_heads.cpp
:language: cpp
```

Next, [step-09](step-09.md): the option surface — the same network emitted several ways, with the
differences in the generated C++ side by side.
