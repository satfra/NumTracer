# Dirac traces

The running example's numerator is a Dirac trace, $T_{\text{num}} = \mathrm{tr}[\slashed p\,\gamma^\mu\,\slashed q\,\gamma^\nu]\,P_{\mu\nu}(l)$.
This page shows the core of it in **two ways**: first, as an explicit 4×4 matrix product (the most concrete view), then the optimized way via the engine's token API (the form the code generator uses).

A gamma trace is a closed loop of 4×4 matrices over spinor indices: the matrices are multiplied in order and the loop is closed by summing the first and last spinor index. A slashed momentum $\slashed p = p^\mu\gamma_\mu$ is one such matrix, carrying the momentum components.
NumTracer multiplies the matrices rather than enumerating the $(2n-1)!!$ index pairings of Wick's theorem, i.e. explicitly using the Clifford algebra to derive a symbolic polynomial in the momentum components.

## 2a — as a matrix product

`Tutorials/02-dirac-traces/trace_raw.cpp` builds each slashed momentum as a 4×4 matrix whose
entries are *polynomials* in the momentum components (`numeric::MPoly`), multiplies, and traces.

```cpp
#include "numtracer/numeric/mpoly.hpp"     // MPoly, Mat4, slashC, matmul, mtrace, eval

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace nm = numtracer::numeric;
using numtracer::Cx;

int main() {
  // The 8 components of two momenta p, q are symbols 0..7: MPoly::var(nsym, i) is "symbol i".
  const int nsym = 8;
  std::array<nm::MPoly, 4> p, q;
  for (int mu = 0; mu < 4; ++mu) {
    p[mu] = nm::MPoly::var(nsym, mu);     // p_mu = symbol mu
    q[mu] = nm::MPoly::var(nsym, 4 + mu); // q_mu = symbol 4+mu
  }

  // slashC(nsym, comp) = sum_mu comp[mu] * gamma^mu : the slashed momentum as a 4x4 of MPoly.
  nm::Mat4 ps = nm::slashC(nsym, p);
  nm::Mat4 qs = nm::slashC(nsym, q);

  // tr(p/ q/) is the trace of the matrix product — a polynomial in the 8 symbols.
  nm::MPoly tr = nm::mtrace(nm::matmul(ps, qs));

  // Evaluate at one point and compare to the closed form 4 p.q.
  std::vector<double> x = {1.0, 0.5, -0.3, 0.2, 0.8, -0.4, 1.2, 0.1};
  Cx got = nm::eval(tr, x, /*no inverse atoms*/ {});

  double pq = 0;
  for (int mu = 0; mu < 4; ++mu) pq += x[mu] * x[4 + mu];

  std::printf("tr(p/ q/)  = %g + %gi   (%zu monomials)\n", got.re, got.im, tr.size());
  std::printf("4 (p.q)    = %g\n", 4.0 * pq);

  const bool ok = std::fabs(got.re - 4.0 * pq) < 1e-12 && std::fabs(got.im) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
```

```bash
cmake --build build --target trace_raw && ./build/trace_raw
```
```text
tr(p/ q/)  = 1.04 + 0i   (4 monomials)
4 (p.q)    = 1.04
ALL TESTS PASSED
```

The pieces:

| symbol                  | meaning                                                                                          |
| ----------------------- | ----------------------------------------------------------------------------------------------- |
| `nsym`                  | the **fixed size of the symbol space** — how many runtime symbols every polynomial is built over (here 8: $p_0..p_3$, $q_0..q_3$) |
| `MPoly::var(nsym, i)`   | the polynomial equal to **symbol `i`** ($x_i$, coefficient 1). `nsym` is the symbol-space size; `i` is the 0-based index of which symbol, so `0 ≤ i < nsym` |
| `slashC(nsym, comp)`    | $\slashed p = \sum_\mu \mathtt{comp}[\mu]\,\gamma^\mu$ as a 4×4 matrix of `MPoly`                |
| `matmul` / `mtrace`     | ordinary 4×4 matrix product and trace, over `MPoly` entries                                      |
| `eval(tr, x, {})`       | plug the symbols' numeric values `x` (one per symbol, so `x.size() == nsym`) into the polynomial (no inverse atoms here) |

```{admonition} nsym is a dimension you fix up front, not a running count
:class: important
The first argument to `MPoly::var(nsym, i)` (and to `MPoly::constant`, `slashC`, …) is **not** "how
many variables I have created so far." It is the fixed dimension of the variable space: internally
every monomial stores an exponent vector of length `nsym`, one slot per symbol. You choose `nsym`
once for the whole computation (here 8), then *address* symbols by index `i` — `var(8, 0)` is
$x_0$, `var(8, 7)` is $x_7$. There is no incremental "add a variable" step; symbol `i` exists for
any `0 ≤ i < nsym` whether or not you use it, and the value of `x_i` is supplied later, positionally,
in the `eval(tr, x, …)` array (`x[i]`).

**Every `MPoly` that is combined (added or multiplied) must share the same `nsym`.** The `+` and `*`
operators assume both operands' exponent vectors have the same length and walk them slot-for-slot;
mixing two different `nsym` values reads exponents out of bounds — undefined behaviour, not a checked
error — and like terms silently fail to combine. So pick `nsym` once, thread the same value through
every builder, and size the `eval` array to match. (When momenta are plain *numbers* rather than
symbols — as in [the full diagram](full-diagram.md) — you set `nsym = 0` and build the components
with `MPoly::constant`, i.e. an empty symbol space.)
```

The trace *is* a matrix product. NumTracer never "knows" the identity
$\mathrm{tr}(\slashed p\,\slashed q) = 4\,p\cdot q$ — it falls out of multiplying the tabulated
gammas, and the surviving polynomial (4 monomials: $4\sum_\mu p_\mu q_\mu$) is exactly $4\,p\cdot q$.

## 2b — via the token API

In practice you don't build matrices by hand: you **describe** the closed chain as a list of
tokens and let `numeric_value` contract it — the form the code generator emits.
`Tutorials/02-dirac-traces/trace_tokens.cpp`:

```cpp
#include "numtracer/numeric/numeric_contract.hpp" // numeric_value, NNet/NTerm, nmet, dslash, dgamma
// … MPoly comp table for p (0..3), q (4..7) as in 3a …

enum { mu, nu }; // name the free gamma legs' Lorentz indices (one unscoped enum, all distinct)

// (1) tr(p/ q/): a closed chain of two slashed momenta.
// dslash({{c, vid}}) = a slashed momentum with components c * (momentum vid).
net::DiracNet chainPQ = {net::dslash({{1.0, 0}}), net::dslash({{1.0, 1}})};
nm::MPoly trPQ = nm::numeric_value(nsym, chainPQ, /*lorentz*/ {}, comp, /*atomDen*/ {});

// (2) tr(gamma^mu p/ gamma_mu q/): two FREE gamma legs (dgamma) on Lorentz indices mu, nu,
// tied together by a metric. A Lorentz network (NNet) is a sum of terms; here one term,
// coefficient 1, with the single factor nmet(mu, nu) = delta_{mu nu}.
net::DiracNet chainG = {net::dgamma(mu), net::dslash({{1.0, 0}}),
                        net::dgamma(nu), net::dslash({{1.0, 1}})};
nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(mu, nu)}}};
nm::MPoly trG = nm::numeric_value(nsym, chainG, lor, comp, {});
```

```bash
cmake --build build --target trace_tokens && ./build/trace_tokens
```
```text
tr(p/ q/)           = 1.04   (= 4 p.q = 1.04, 4 monomials)
tr(g^mu p/ g_mu q/) = -2.08   (= -2 tr(p/ q/) = -2.08)
ALL TESTS PASSED
```

The tokens:

| token                | meaning                                                                                                                            |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `dslash({{c, vid}})` | a slashed momentum whose components are $c\times$ momentum `vid` (a linear combination, so a propagator $l-p$ is `{{1,l},{-1,p}}`) |
| `dgamma(id)`         | a free $\gamma^{id}$ with an **open** Lorentz index, to be contracted by a Lorentz network                                         |
| `nmet(a, b)`         | the metric $\delta_{ab}$, one factor of a Lorentz network term                                                                     |
| `NNet` / `NTerm`     | a Lorentz network: a sum of terms, each a coefficient times a product of factors                                                   |

Chain (1) has no free legs, so it is a scalar — the same $4\,p\cdot q$ as 3a. Chain (2) has two
free gamma legs contracted by the metric, giving $\mathrm{tr}(\gamma^\mu\slashed p\,\gamma_\mu
\slashed q) = -2\,\mathrm{tr}(\slashed p\,\slashed q) = -8\,p\cdot q$ (the 4D identity
$\gamma^\mu\slashed a\,\gamma_\mu = -2\slashed a$). The `MPoly` result is symbolic — `numeric_value`
returns the *polynomial*, and `eval` plugs numbers in afterward.

Next, [Lorentz networks](lorentz-networks.md): the projector $P_{\mu\nu}(l)$ that contracts those
free legs.
