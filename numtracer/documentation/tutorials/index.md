# Tutorials

These tutorials follow one running example — the **quark self-energy**, the one-loop diagram in
which a quark emits and reabsorbs a gluon — starting with the tensor network itself, all the way to a generated,
validated kernel that contracts the network.

Each page leads with a complete, runnable program. You can write and compile them yourself, or find and compile them in the repository's `Tutorials/` folder. The programs depend only on `NumTracer` itself and require a C++20 compiler.

## The running example

Deriving the flow of the quark two-point function (two functional derivatives of the Wetterich equation w.r.t. quarks) produces exactly **two diagrams**: the same quark–gluon triangle, once with the regulator insertion $\partial_t R$ on the **gluon** line and once on the **quark** line.
Projecting onto the wave function and performing the colour and Dirac algebra, the integrand factorises into

$$
\text{integrand} \;=\; \underbrace{C_F}_{\text{colour}} \times
\Big(\,\text{term}_1\big|_{q=l-p} \;+\; \text{term}_2\big|_{q=l+p}\Big),
$$

where both diagrams carry the same colour structure $T^aT^a = C_F = 4/3$ for SU(3), and each term carries the Dirac-trace numerator

$$
T_{\text{num}} \;=\; \mathrm{tr}\!\big[\slashed p\,\gamma^\mu\,\slashed q\,\gamma^\nu\big]\,
P_{\mu\nu}(l),
\qquad
P_{\mu\nu}(l) = \delta_{\mu\nu} - \frac{l_\mu l_\nu}{l^2}.
$$

The colour factor and the trace $T_{\text{num}}$ are the two pieces NumTracer computes; the rest (dressings, regulator, etc.) is provided by the user in ordinary C++.

```{admonition} Prerequisite
:class: tip
Read [**Key concepts**](../getting_started/concepts.md) first — it defines how a network is described (indices contract when they share a label; the frame fixes each momentum's components; scalar products are the kernel's symbols).
```

## Building the programs

To directly compile the programs presented in the tutorials, you can find them in the repository's `Tutorials/` folder. There is a standalone CMake project that consumes NumTracer like any external user would:

```bash
cd Tutorials
cmake .                # configure once
make -j4               # build every tutorial
ctest --test-dir build # run them all as a self-checking test suite
```

Each page below also gives the one-target command (`make <name>`) and the program's real output.

## The path through these pages

1. **[Colour factors](color-contraction.md)** — folding a momentum-free colour network to a single number, $C_F = 4/3$.
2. **[Dirac traces](dirac-traces.md)** — the trace $T_{\text{num}}$ hand-coded two ways: first as an explicit 4×4 matrix product, then via the engine's token API.
3. **[Lorentz networks](lorentz-networks.md)** — contract a network of metrics, vectors and projectors to a scalar polynomial over the loop frame.
4. **[A full diagram](full-diagram.md)** — assemble the whole $T_{\text{num}}$ by hand and validate it two ways (numeric engine, closed form).
5. **[Generating kernels from Mathematica](generating-kernels.md)** — the front-end that writes all of the above from a tiny DSL.
6. **[Finite temperature](finite-temperature.md)** — the running example in the heat-bath frame: electric/magnetic projectors, the spatial scalar product, and a distinct Matsubara component.
7. **[Per-flavour and per-component dressings](dressed-flavour.md)** — dress individual flavours or colour/flavour components differently: the per-flavour split and the group-diagonal dressing heads.

For the design behind the engine, see [Under the Hood](../internals/index.md).

```{toctree}
:maxdepth: 1

color-contraction
dirac-traces
lorentz-networks
full-diagram
generating-kernels
finite-temperature
dressed-flavour
```
