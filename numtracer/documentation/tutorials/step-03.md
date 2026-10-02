# step-03: Dirac traces — matrices, then tokens

*Builds on: [step-01](step-01.md) · Built on by: [step-05](step-05.md), [step-18](step-18.md) ·
Tags: `dirac`, `trace` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

The Dirac sector is where the naive approach stops being merely wasteful and becomes impossible, so
it is worth being precise about what a gamma trace *is* before writing any code.

### What a trace actually is

A quark propagating through a diagram carries a spinor index. Every vertex it passes inserts a
$4\times4$ matrix on that index; when the quark line closes into a loop, the first and last spinor
indices are identified and summed. That sum-over-the-diagonal is the trace. So

$$
\mathrm{tr}\big[\slashed p\,\slashed q\big]
$$

is *literally* the trace of a product of two $4\times4$ matrices, where a **slashed momentum**
$\slashed p = p^\mu\gamma_\mu$ is the matrix you get by weighting the four gammas with the four
components of $p$.

The textbook way to evaluate such a trace is Wick's theorem: sum over all pairings of the gamma
indices, with a sign per crossing. For $2n$ gammas that is $(2n-1)!!$ terms — 3 for four gammas,
10395 for fourteen, and a four-point vertex flow reaches well past that. Worse, the intermediate is
a sum of scalar-product monomials that must then be collected.

NumTracer does not use Wick's theorem. It **multiplies the matrices**. The gammas are typed out as
`constexpr` tables in `dirac/dirac_data.hpp` (Hermitian, chiral/Weyl basis, Euclidean, so
$\{\gamma^\mu,\gamma^\nu\} = 2\delta^{\mu\nu}$), each with only four nonzero entries out of
sixteen. Multiplying a chain is $O(n)$ matrix products of a fixed $4\times4$ size, and the Clifford
algebra is *implied* rather than applied — every identity you know falls out of the arithmetic.

```{admonition} The entries are polynomials, not numbers
:class: important
This is the trick that makes the matrix-product approach work symbolically. The engine does not
multiply matrices of `double`; it multiplies matrices of `nt::Poly` — multivariate polynomials in the
frame's symbols. So $\slashed p$ is a $4\times4$ array whose entries are things like
"$p_0 + i p_3$", and the product of two such matrices has entries that are polynomials of degree 2.
The trace is then one polynomial, already collected, with no separate simplification pass.
```

This step shows the same trace twice: first with the matrices explicit, so you can see there is no
magic; then through the **token API**, which is how every generated kernel expresses it.

## The commented program, part a — as a matrix product

`Tutorials/step-03-dirac-traces/trace_raw.cpp`.

### The symbol space

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_raw.cpp
:language: cpp
:start-after: "@snip begin: symbols"
:end-before: "@snip end: symbols"
```

The eight components are the frame's symbols, $p_\mu$ and $q_\mu$. No frame has been chosen
— the eight are independent — which is the right setting to check an algebraic identity: if the
polynomial identity holds for independent symbols, it holds in every frame. A `Symbol` converts to
the polynomial it stands for, and symbols combine with `+`, `-`, `*` and numbers.

```{admonition} Symbols are declared first, then frozen
:class: important
Every polynomial of a frame stores one exponent per frame symbol, so the symbol list must be
complete before the first polynomial exists. The frame therefore *freezes* its symbol list the
first time a polynomial is made from it (a momentum, an arithmetic expression, `slashC`). Declaring
a symbol after that throws. Values are supplied at evaluation time, in declaration order:
`F.at(p0, p1, …)`, or by name: `F.at({{P0, 1.0}, …})`.
```

### Slash, multiply, trace

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_raw.cpp
:language: cpp
:start-after: "@snip begin: slash"
:end-before: "@snip end: slash"
```

Three lines, and every one is ordinary linear algebra. `F.slashC(k)` builds
$\sum_\mu k_\mu\gamma^\mu$ as a `Mat4` of polynomials; `matmul` is the $4\times4$ product;
`mtrace` sums the diagonal. There is no Dirac-algebra code path here at all — `matmul` does not
know it is multiplying gammas.

## Results, part a

```bash
cmake --build build --target trace_raw && ./build/trace_raw
```

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_raw.expected.txt
:language: text
```

Two things to notice.

**The identity was never coded.** NumTracer contains no rule saying
$\mathrm{tr}(\slashed p\slashed q) = 4\,p\cdot q$. It multiplied two matrices of polynomials and
traced. The result agrees because the tabulated gammas satisfy the Clifford algebra — which is
itself checked, in `tests/`, against a runtime oracle, so a typo in the tables fails a test rather
than silently producing wrong physics.

**"4 monomials" is the whole answer.** The polynomial is $4\sum_\mu p_\mu q_\mu$ — four terms, one
per $\mu$, each with coefficient 4. That *is* $4\,p\cdot q$, written in components. And the
imaginary part is exactly zero, not $10^{-17}$: the individual matrix entries are complex (the
gammas have $\pm i$ in them) and the imaginary parts cancel term by term in exact arithmetic.

## The commented program, part b — via the token API

In practice you never build the matrices. You **describe** the closed chain as a list of tokens and
let the engine contract it — which is exactly the form the code generator emits.

```{admonition} What `F.trace` does
:class: note
`F.trace(chain, net)` *contracts*: it multiplies and sums over every shared index.

1. It closes the `chain` into $4\times4$ gamma products and takes the spinor trace (part a's
   `matmul`/`mtrace`, done for you). Legs left open with `gamma(mu)` survive as *free Lorentz
   indices* on the resulting tensor, and $\mathrm{tr}\,\mathbb{1} = 4$ rides along.
2. It contracts that tensor against the Lorentz network `net`, summing every index shared between a
   `gamma` leg and a network factor.
3. It returns the surviving scalar as one polynomial.

`net` may be left out when the chain has no open legs. Every index must end up contracted; an
index that occurs once is an error.
```

`Tutorials/step-03-dirac-traces/trace_tokens.cpp`. The frame has the same eight symbols, now as the components of two momenta:

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.cpp
:language: cpp
:start-after: "@snip begin: comp"
:end-before: "@snip end: comp"
```

### A chain with no free legs

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.cpp
:language: cpp
:start-after: "@snip begin: chain-pq"
:end-before: "@snip end: chain-pq"
```

The chain is a `DiracChain`: a list of tokens **in trace order**, implicitly cyclic — the last
spinor index is tied back to the first, because that is what closing a quark loop means. So this
list of two `slash` tokens already *is* $\mathrm{tr}(\slashed p\slashed q)$, and no Lorentz
network is needed because nothing is left open.

`nt::slash(p)` is $\slashed p$. Momenta combine like vectors, so an internal propagator momentum is
written `nt::slash(l - p)` — no new momentum is declared for it; the frame knows only the
independent ones.

### A chain with free legs

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.cpp
:language: cpp
:start-after: "@snip begin: chain-g"
:end-before: "@snip end: chain-g"
```

`nt::gamma(mu)` is a $\gamma^\mu$ whose Lorentz index is **open**. A chain containing open legs is not
a scalar — it is a tensor with those free indices — so it must be handed a Lorentz network that
ties them off. Here that network is a single metric, which contracts $\mu$ with $\nu$ and gives
$\mathrm{tr}(\gamma^\mu\slashed p\,\gamma_\mu\slashed q)$.

This is the shape of every real diagram: the Dirac chain carries the gamma structure with open
gluon legs, and the Lorentz network carries the gluon propagators and projectors that close them.
[step-05](step-05.md) assembles exactly that.

## Results, part b

```bash
cmake --build build --target trace_tokens && ./build/trace_tokens
```

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.expected.txt
:language: text
```

The first line reproduces part a exactly, as it must — the token path and the matrix path are the
same computation with different bookkeeping.

The second is the 4-dimensional identity $\gamma^\mu\slashed a\,\gamma_\mu = -2\slashed a$, giving
$-2\,\mathrm{tr}(\slashed p\slashed q) = -8\,p\cdot q$. The factor $-2$ is where the Euclidean,
strictly-4-dimensional convention shows: in $d$ dimensions it would be $(2-d)$, and there is no
general-$d$ mode in NumTracer. If your calculation needs dimensional regularisation, this is the
boundary — see [Scope & conventions](../getting_started/scope-and-conventions.md).

Note finally that `F.trace` returned a *polynomial*, and `F.eval` put numbers in afterwards.
That separation is the whole basis of code generation: the polynomial is computed once at build
time, and the kernel that ships evaluates it millions of times.

## Possibilities for extensions

1. **Odd numbers of gammas.** Trace a chain of three slashes. You should get exactly zero — every
   trace of an odd number of gammas vanishes. Confirm the polynomial has *no* monomials at all
   (`tr.size() == 0`), not merely monomials that evaluate to zero.

2. **The four-gamma identity.** Verify
   $\mathrm{tr}(\slashed a\slashed b\slashed c\slashed d) = 4[(a\!\cdot\!b)(c\!\cdot\!d) -
   (a\!\cdot\!c)(b\!\cdot\!d) + (a\!\cdot\!d)(b\!\cdot\!c)]$ with four independent momenta
   (sixteen symbols). Count the monomials before you run it and see whether your guess was right.

3. **$\gamma_5$.** Add `nt::gamma5()` to a chain of four slashes. The result is proportional to the
   Levi-Civita tensor and vanishes unless all four momenta are linearly independent — so it is
   *zero* in the 1-angle frames of [step-04](step-04.md) and nonzero with four generic momenta.
   Check both.

4. **Break it deliberately — forget to close the legs.** Trace the chain of `trG` without its
   metric. An open Lorentz index with nothing to contract it against is a real error, and the
   engine says so (the message is quoted in full by `tests/test_frame.cpp`):

   ```text
   numtracer: Lorentz index id 0 is OPEN (occurs once) — ...
   ```

   The engine cannot catch the other mistake: a metric on the *wrong* pair of labels, as long as
   every label still occurs twice, is a different but perfectly valid contraction.

5. **Measure the scaling.** Time chains of 4, 8, 12, 16 slashes. The matrix-product cost grows
   linearly in the chain length; the Wick pairing count grows as $(2n-1)!!$. Plot both on a log
   axis and you have the argument for the whole design.

## The plain programs

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_raw.cpp
:language: cpp
```

```{literalinclude} ../../../Tutorials/step-03-dirac-traces/trace_tokens.cpp
:language: cpp
```

Next, [step-04](step-04.md): the Lorentz network that ties those open legs off, and what happens to
the $1/l^2$ it drags along.
