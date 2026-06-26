# Sector data

> Headers: `dirac/dirac_data.hpp`, `dirac/dense_trace.hpp`, `sun/sun_data.hpp`, `sun/sun.hpp`,
> `dense/dtensor.hpp`

The contraction engine is ignorant of physics — it folds matrices and sums indices. The physics
constants live in small, typed-out tables that the engine reads. Each sector has a `constexpr`
**data table** of its constants, cross-checked against a runtime **oracle** so a typo fails a
test rather than silently giving wrong physics.

## Dirac: gamma matrices

The Dirac sector handles the spin-½ algebra: 4×4 matrices acting on spinor indices. The
generators are the gamma matrices $\gamma^\mu$ ($\mu = 1..4$ in 4D Euclidean space), satisfying
the Clifford algebra $\{\gamma^\mu, \gamma^\nu\} = 2\,\delta^{\mu\nu} I$. NumTracer uses the
Hermitian gammas in the chiral (Weyl) basis:

$$
\gamma^k = \begin{pmatrix} 0 & -i\sigma^k \\ i\sigma^k & 0 \end{pmatrix}\ (k=1,2,3),
\qquad
\gamma^4 = \begin{pmatrix} 0 & I \\ I & 0 \end{pmatrix},
\qquad
\gamma_5 = \mathrm{diag}(1,1,-1,-1).
$$

`dirac/dirac_data.hpp` types these out as `constexpr` tables of `Cx`, indexed
`kGamma[mu][i][j]` $= (\gamma^\mu)_{ij}$. Each gamma has only four nonzero entries out of
sixteen; that sparsity, and the block-antidiagonal Weyl structure, are what the trace exploits.

### The chiral dense trace

`dirac/dense_trace.hpp` provides a sparsity-aware dense gamma trace, `chiral_gamma_trace`. A
slashed momentum in the Weyl basis is block-antidiagonal, $\slashed p = \big[\begin{smallmatrix}0
& P\\ Q & 0\end{smallmatrix}\big]$, so a chain of them folds into products of 2×2 chiral blocks
rather than 4×4 matrices — a few times faster than the brute 4×4 chain — and yields
$\mathrm{tr}(\text{odd}) = 0$ structurally. This is the same idea the
[numeric engine](numeric-engine.md) uses when it traces a closed Dirac chain by matrix products.

## SU(N): colour algebra

The colour sector implements SU(N) (N = 3 for QCD). Its objects are the generators $T^a$ — the
$N\times N$ Hermitian traceless matrices of the fundamental representation, $a = 1..N^2-1$ — and
the structure constants $f^{abc}$ defined by $[T^a, T^b] = i f^{abc} T^c$. Colour factors are
contractions of these, e.g. $f^{acd} f^{bcd} = N\,\delta^{ab}$ and
$\mathrm{tr}(T^a T^b) = \tfrac12 \delta^{ab}$.

`sun/sun_data.hpp` types out the SU(2)/SU(3) tables, and `sun/sun.hpp` provides the runtime
**oracle** `SUNBuilder<N>`: its constructor builds the generators (generalized Gell-Mann construction)
and the structure constants $f^{abc} = -2i\,\mathrm{tr}([T^a,T^b]T^c)$ for any `N`, which is the
source of truth the typed-out tables are checked against. A colour network with no free indices
folds to a single number — the [numeric engine](numeric-engine.md) does this at build time over
these tables (`network/sun_net.hpp`), so the kernel never carries a colour tensor. A network that
dresses each colour/flavour component differently (a group-diagonal $\delta$) folds instead with
`sun_value_dressed` to a polynomial $\sum_a c_a Z_a$ over runtime per-component dressing leaves;
`sun_value_cx` is left untouched, so colour-blind flows are unchanged.

`sun.hpp` also defines `Mat<N>`, `matmul`, and `trace` — a plain dense complex matrix type used
by the oracle and the tests.

## The dense numeric tensor

`dense/dtensor.hpp` defines `dense::DTensor`, a dense numeric tensor over labelled axes
`Ax<Id, Dim>`. It contracts any network — Lorentz, Dirac, SU(N), or a mix — entry-for-entry with
runtime `Cx` numbers and **no** structural pruning. That makes it the entry-for-entry **oracle**
the generated kernels are validated against, and the target of the codegen's `"Dense"` backend
(see [codegen](codegen.md)). It and the symbolic contraction share the same `Ax<Id, Dim>` axis
language and a division-free index *odometer* (inline in `dense/dtensor.hpp`), a mixed-radix counter
that walks every index combination while maintaining running flat offsets — so the inner contraction
loop needs no division or modulo and vectorises well.
