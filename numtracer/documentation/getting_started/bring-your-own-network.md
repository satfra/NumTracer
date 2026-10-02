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
```text
a . b        = 0.26   (expect 0.26)
a . P(k) . a = 0.403649   (expect a^2 - (a.k)^2/k^2 = 0.403649)
ALL TESTS PASSED
```

That is the whole engine in miniature: describe a network as a list of factors sharing integer
labels, hand it a frame that fixes each vector's components, and read off the scalar. Add a Dirac
chain or SU($N$) factors and nothing else changes — the sectors compose (see the
[tutorials](../tutorials/index.md)).

## Two ways to drive it

| Path | You write | Needs | Good for |
|---|---|---|---|
| **C++ API** | the network with `nvec`/`nmet`/`nprojT`/`dgamma`/`SUN::T`… and call `env.numeric_value` on a `LorentzEnv` | nothing but a C++20 compiler — no external libraries | hand-built traces, embedding in your own code, the numeric oracle |
| **Mathematica DSL** | the network with `ntVec`/`ntMetric`/`ntTransProj`/… and call `NumTrace` + `MakeNTKernel` | a Wolfram kernel | generating a committed, lowered C++ kernel from a symbolic network |

Two things are worth stating plainly, because it is easy to assume otherwise:

- **What FunKit is needed for.** `FromFunKit` is a convenience importer for flows *already derived*
  in the FunKit/DiFfRG toolchain; you do not need it — hand-build the DSL network and give it to
  `NumTrace` directly, as [step-06](../tutorials/step-06.md) does. `MakeNTKernel`, however, writes the
  kernel through FunKit's C++ emitter, so the Mathematica code generator needs FunKit installed
  either way. The C++ API needs neither.
- **FORM is not needed** to use NumTracer or to generate your own kernels. It appears only when
  *regenerating the project's own reference-test oracles*, never on your path.

## The dictionary: DSL head ↔ C++ builder

The two front-ends describe the *same* engine, so every DSL head has a C++ builder:

| Sector | Mathematica DSL | C++ builder |
|---|---|---|
| metric $\delta_{\mu\nu}$ | `ntMetric[mu, nu]` | `nm::nmet(mu, nu)` |
| vector $q_\mu$ | `ntVec[q, mu]` | `nm::nvec(mu, {{coeff, vid}, …})` |
| transverse projector | `ntTransProj[q, mu, nu]` | `nm::nprojT(mu, nu, {{1.0, vid}}, atom)` |
| longitudinal projector | `ntLongProj[q, mu, nu]` | `nm::nprojL(mu, nu, {{1.0, vid}}, atom)` |
| Levi-Civita $\varepsilon$ | (γ5 trace) | `nm::neps(a, b, c, d)` |
| gamma $\gamma^\mu$ (open leg) | `ntGamma[mu, …]` | `net::dgamma(mu)` |
| slashed $\slashed q$ | slashed momentum | `net::dslash({{coeff, vid}, …})` |
| SU($N$) generator $(T^a)_{ij}$ | `ntSUNT[N, a, i, j]` | `SUN::T(N, a, i, j)` |
| SU($N$) structure constant | `ntSUNf[N, a, b, c]` | `SUN::f(N, a, b, c)` |
| adjoint / fund. $\delta$ | `ntSUNDeltaAdj` / `ntSUNDeltaFund` | `SUN::deltaAdj` / `SUN::deltaFund` |

(The `vid` is a vector id into your component table; the `atom` is the id of a projector's
$1/q^2$. Colour networks are built as a `SUNNet` and folded with `network::sun_value_cx` — see the
[step-02](../tutorials/step-02.md).)

## Mapping your own theory

1. **List your indices and give each a label.** Use one unscoped `enum` across *all* sectors so
   every label is distinct by construction (`enum { mu, nu, a, b, A, B };`). Reusing a label is how
   you say "sum these together"; a label used once stays free.
2. **Pick a frame.** Choose reference components for each vector — for a one-angle loop, one vector
   along an axis and another at an angle (see [Key concepts](concepts.md#runtime-numbers-come-from-the-frame)).
   Keep components you want as runtime symbols symbolic (`env.var`, from a `LorentzEnv env(nsym)`);
   the rest are `env.constant`.
3. **Build the network and contract.** Assemble the Dirac chain (`DiracNet`), the Lorentz network
   (`NNet`), and any colour network, then call `env.numeric_value` (Lorentz+Dirac) and `sun_value_cx`
   (colour). You get back a scalar `MPoly` in your frame symbols.
4. **Read it or lower it.** `eval` the polynomial at a point, or hand it to the codegen to emit a
   flat C++ kernel ([step-05](../tutorials/step-05.md)).

If any object or convention above does not match your problem, check
[Scope & conventions](scope-and-conventions.md) — that page is the exact boundary of what the
engine represents.
