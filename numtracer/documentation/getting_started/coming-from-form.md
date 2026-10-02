# Coming from FORM

If you have done Dirac traces in FORM, most of NumTracer will look familiar: you declare vectors and
indices, write a gamma chain, ask for its trace. This page maps the vocabulary in one table, shows
the same trace three ways, and lists what is genuinely different.

## The same trace, three ways

$$
\mathrm{tr}\big(\gamma^\mu\slashed p\,\gamma_\mu\slashed q\big) \;=\; -8\,p\cdot q
$$

**FORM:**

```text
Vectors p, q;
Indices mu;
Local T = g_(1, mu, p, mu, q);
trace4, 1;
Print;
.end
```

FORM answers with `T = -8*p.q;` — a polynomial in the scalar product `p.q`.

**NumTracer, C++** (`namespace nt = numtracer;`):

```cpp
nt::Frame F;                                      // the momenta need components
auto P0 = F.symbol("p0"), P1 = F.symbol("p1"), P2 = F.symbol("p2"), P3 = F.symbol("p3");
auto Q0 = F.symbol("q0"), Q1 = F.symbol("q1"), Q2 = F.symbol("q2"), Q3 = F.symbol("q3");
nt::Momentum p = F.momentum(P0, P1, P2, P3);
nt::Momentum q = F.momentum(Q0, Q1, Q2, Q3);
auto [mu, nu] = F.indices<2>();

nt::Poly T = F.trace({nt::gamma(mu), nt::slash(p), nt::gamma(nu), nt::slash(q)}, nt::metric(mu, nu));
double v = F.eval(T, F.at(1.0, 0.5, -0.3, 0.2, 0.8, -0.4, 1.2, 0.1)).re;   // = -8 p.q at this point
```

NumTracer answers with a polynomial in the momentum *components*,
$-8(p_0q_0 + p_1q_1 + p_2q_2 + p_3q_3)$, which `F.eval` turns into a number.
([step-03](../tutorials/step-03.md) is this program, run.)

**NumTracer, Mathematica DSL:**

```mathematica
net = ntGamma[mu, d1, d2] ntVec[p, rho] ntGamma[rho, d2, d3] ntGamma[mu, d3, d4] ntVec[q, sig] ntGamma[sig, d4, d1];
ntk = NumTrace[net, "Frame" -> frame, "Args" -> args];   (* then MakeNTKernel[ntk, …] writes C++ *)
```

Spinor indices `d1 → d2 → d3 → d4 → d1` close the loop; a slash is a gamma contracted with a vector.
([step-06](../tutorials/step-06.md) is a complete script.)

## The vocabulary

| FORM | NumTracer C++ | NumTracer DSL |
|---|---|---|
| `Vectors p;` | `nt::Momentum p = F.momentum(c0, c1, c2, c3);` — components required | a momentum symbol, placed by the frame |
| `Symbols x;` | `nt::Symbol X = F.symbol("x");` | any symbol |
| `Indices mu, nu;` | `auto [mu, nu] = F.indices<2>();` | any symbol used as a label |
| `p(mu)` | `nt::vec(mu, p)` | `ntVec[p, mu]` |
| `d_(mu, nu)` | `nt::metric(mu, nu)` | `ntMetric[mu, nu]` |
| `e_(mu, nu, rho, si)` | `nt::epsilon(mu, nu, rho, si)` | (produced by $\gamma_5$ traces) |
| `p.q` | `F.contract(nt::vec(mu, p) * nt::vec(mu, q))` | `ntSP[p, q]` |
| `g_(1, mu)` | `nt::gamma(mu)` in a `nt::DiracChain` | `ntGamma[mu, d1, d2]` |
| `g_(1, p)` | `nt::slash(p)` | `ntVec[p, mu] ntGamma[mu, d1, d2]` |
| `g5_(1)` | `nt::gamma5()` | `ntGamma5[d1, d2]` |
| `trace4, 1;` | `F.trace(chain, net)` | `NumTrace` (then `MakeNTKernel`) |
| `Local T = …;` | `nt::Poly T = …;` | `net = …;` |
| `id …;` (substitutions) | — none: the frame fixes every component | — none |
| `Print;` | `F.eval(T, F.at(…))` at a point | the generated kernel |
| colour (`color.h`) `T(i, j, a)`, `f(a, b, c)` | `nt::SUN su3(3);` `su3.T(a, i, j)`, `su3.f(a, b, c)` | `ntSUNT[3, a, i, j]`, `ntSUNf[3, a, b, c]` |

The transverse projector $\delta_{\mu\nu} - k_\mu k_\nu/k^2$, which FORM users write out by hand, is
one object here: `nt::projT(mu, nu, k)`.

## What is different

**The answer is numeric in the kinematics, not symbolic in scalar products.** FORM manipulates
`p.q` as a symbol. NumTracer gives every momentum explicit components in a *frame* and contracts
over them, so a scalar product never appears: the result is a polynomial in the frame's symbols
(magnitudes, angle cosines). The payoff is speed and a direct path to a C++ kernel; the price is
that you choose a frame. For a loop integrand the natural frame is cheap and exact
([step-04](../tutorials/step-04.md)).

**Four dimensions, Euclidean.** Index sums run over 0–3 and the metric is $\delta_{\mu\nu}$, so
$\gamma^\mu\gamma_\mu = 4$ and $\gamma^\mu\slashed a\gamma_\mu = -2\slashed a$. There is no
`tracen` and no dimensional regularisation. See [Scope & conventions](scope-and-conventions.md).

**No pattern matching.** There are no `id` statements and no user-defined functions acting on
expressions; the network you write is the whole input, and the algebra is fixed. Things FORM users
do with `id` — substituting a propagator, a projector, a basis element — are written into the
network directly.

**The output is a kernel.** A FORM trace ends in a printed expression; a NumTracer trace ends in a
polynomial that is lowered to straight-line C++ ([step-05](../tutorials/step-05.md)), or, through
the Mathematica front-end, in a complete kernel for an integrator ([step-06](../tutorials/step-06.md)).

**SU($N$) is exact and separate.** An SU($N$) factor (colour, flavour, …) folds to an exact number,
computed apart from the Dirac/Lorentz trace and multiplied in ([step-02](../tutorials/step-02.md)).
