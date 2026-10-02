# Bring your own network

The worked examples in this guide are QCD/fRG loop integrands, but that is just what the authors
needed traced — the engine itself is theory-agnostic. This page is the general on-ramp: a complete
program with **no physics vocabulary**, the two ways to drive the engine, and the dictionary that
maps *your* tensors onto its heads. If you have a network of metrics, vectors, projectors, gamma
matrices, and/or SU($N$) factors — in the [conventions](scope-and-conventions.md) NumTracer
assumes — this is how you contract it.

## Hello, tensor network

The smallest possible use, in pure C++ against the library. It contracts two index
networks over plain 4-vectors and checks each against its closed form:

$$
a\cdot b = a_\mu\,\delta^{\mu\nu}\,b_\nu,
\qquad
a\cdot P(k)\cdot a = a^2 - \frac{(a\cdot k)^2}{k^2},
\quad P(k)_{\mu\nu} = \delta_{\mu\nu} - \frac{k_\mu k_\nu}{k^2}.
$$

The program is `Tutorials/step-01-hello-network/hello_network.cpp`, walked through line by line in
[**step-01**](../tutorials/step-01.md). Its core is these two networks:

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:start-after: "@snip begin: dot"
:end-before: "@snip end: dot"
```

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:start-after: "@snip begin: proj"
:end-before: "@snip end: proj"
```

```bash
cmake -S Tutorials -B Tutorials/build && cmake --build Tutorials/build --target hello_network
./Tutorials/build/hello_network
```
```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.expected.txt
:language: text
```

That is the whole engine in miniature: declare a frame that fixes each vector's components,
describe a network as a product of factors sharing index labels, and read off the scalar. Add a Dirac
chain or SU($N$) factors and nothing else changes — the sectors compose (see the
[tutorials](../tutorials/index.md)).

## Two ways to drive it

| Path | You write | Needs | Good for |
|---|---|---|---|
| **C++ API** | an `nt::Frame`, the network with `vec`/`metric`/`projT`/`gamma`/`slash`…, and call `F.trace` / `F.contract`; SU($N$) factors through an `nt::SUN` | nothing but a C++20 compiler — no external libraries | hand-built traces, embedding in your own code, checking a result; `to_genprog` + `emit_cpp` lower one polynomial to a C++ function |
| **Mathematica DSL** | the network with `ntVec`/`ntMetric`/`ntTransProj`/… and call `NumTrace` + `MakeNTKernel` | a Wolfram kernel and FunKit | a complete generated kernel: every diagram, the `fill()`, the signature, dressings and SU($N$) factors |

Two things are worth stating plainly, because it is easy to assume otherwise:

- **What FunKit is needed for.** `FromFunKit` is a convenience importer for flows *already derived*
  in the FunKit/DiFfRG toolchain; you do not need it — hand-build the DSL network and give it to
  `NumTrace` directly, as [step-06](../tutorials/step-06.md) does. `MakeNTKernel`, however, writes the
  kernel through FunKit's C++ emitter, so the Mathematica code generator needs FunKit installed
  either way. The C++ API needs neither.
- **FORM is not needed** to use NumTracer or to generate your own kernels. It appears only when
  *regenerating the project's own reference-test oracles*, never on your path.

## The dictionary: DSL head ↔ C++ builder

The two front-ends describe the *same* engine, so every DSL head has a C++ builder (`nt` is
`namespace nt = numtracer;`; `F` a `nt::Frame`, `su` an `nt::SUN`):

| Object | Mathematica DSL | C++ |
|---|---|---|
| momentum $q$ | a symbol placed by the frame, e.g. `propFrame[…]` | `nt::Momentum q = F.momentum(c0, c1, c2, c3);` |
| Lorentz index labels | any symbols | `auto [mu, nu] = F.indices<2>();` |
| metric $\delta_{\mu\nu}$ | `ntMetric[mu, nu]` | `nt::metric(mu, nu)` |
| vector $q_\mu$ | `ntVec[q, mu]` | `nt::vec(mu, q)` |
| transverse projector | `ntTransProj[q, mu, nu]` | `nt::projT(mu, nu, q)` |
| longitudinal projector | `ntLongProj[q, mu, nu]` | `nt::projL(mu, nu, q)` |
| finite-$T$ electric / magnetic projector | `ntElectricProj` / `ntMagneticProj` | `nt::projE` / `nt::projM` |
| Levi-Civita $\varepsilon_{\mu\nu\rho\sigma}$ | (from a $\gamma_5$ trace) | `nt::epsilon(mu, nu, rho, sigma)` |
| product / sum of tensors | `*` / `+` | `*` / `+` |
| $\gamma^\mu$ (open leg) | `ntGamma[mu, d1, d2]` | `nt::gamma(mu)` in a chain |
| slashed $\slashed q$ | `ntVec[q, mu] ntGamma[mu, d1, d2]` | `nt::slash(q)` in a chain |
| $\gamma_5$ | `ntGamma5[d1, d2]` | `nt::gamma5()` in a chain |
| closed spinor loop | spinor labels `d1 → d2 → … → d1` | `nt::DiracChain{…}` in trace order |
| SU($N$) index labels | any symbols | `auto [a] = su.adjoint<1>(); auto [i, j] = su.fundamental<2>();` |
| generator $(T^a)_{ij}$ | `ntSUNT[N, a, i, j]` | `su.T(a, i, j)` |
| structure constant $f^{abc}$ | `ntSUNf[N, a, b, c]` | `su.f(a, b, c)` |
| adjoint / fundamental $\delta$ | `ntSUNDeltaAdj` / `ntSUNDeltaFund` | `su.delta(a, b)` / `su.delta(i, j)` |
| contract | `NumTrace` + `MakeNTKernel` | `F.trace(chain, net)`, `F.contract(net)`, `su.value(net)` |

## Mapping your own theory

1. **List your indices.** Get each sector's labels from its owner: Lorentz labels from the frame,
   SU($N$) labels from one `nt::SUN` object per group. Reusing a label is how you say "sum these
   together"; a label used once is an error.
2. **Pick a frame.** Declare the runtime scalars as symbols and the momenta from them — for a
   one-angle loop, one vector along an axis and another at an angle (see
   [Key concepts](concepts.md#runtime-numbers-come-from-the-frame)). Components that are fixed
   numbers are just numbers.
3. **Build the network and contract.** Assemble the Dirac chain and the Lorentz network and call
   `F.trace(chain, net)` (or `F.contract(net)` without a chain); fold each SU($N$) group with
   `su.value(...)`. You get back a polynomial in your frame's symbols, and the SU($N$) numbers that
   multiply it.
4. **Read it or lower it.** `F.eval(poly, F.at(...))` evaluates the polynomial at a point;
   `nt::to_genprog` + `nt::emit_cpp` lower it to a straight-line C++ function, and `F.emit_fill`
   prints the function that computes its inputs ([step-05](../tutorials/step-05.md)). A complete
   kernel — many diagrams, dressings, the integrator-facing signature — is the job of the
   Mathematica front-end ([step-06](../tutorials/step-06.md)).

If any object or convention above does not match your problem, check
[Scope & conventions](scope-and-conventions.md) — that page is the exact boundary of what the
engine represents.
