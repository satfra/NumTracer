# Colour factors

A colour factor in QCD arises through a contraction of SU(N) generators $T^a$ and structure constants $f^{abc}$ with no free indices and no momenta — so it is usually **just a number**, known before any grid point is evaluated. Both quark-self-energy diagrams exchange one gluon, so both carry the same colour structure

$$
T^a_{ij}\,T^a_{jk} = C_F\,\delta_{ik},
\qquad C_F = \frac{N^2-1}{2N} = \frac{4}{3}\ \text{for SU(3)},
$$

which factors out of the loop.

## Folding a colour network sparsely

The obvious way to evaluate $T^a_{ij}T^a_{jk}$ is to build both generators as full
$(N^2-1)\times N\times N$ tensors of numbers and sum **every entry** of the contraction. That
works, but it does exactly what a colour factor never needs — it *materialises* the whole colour
tensor only to sum it back down to one scalar. It does not scale, either: the colour network of a
four-gluon vertex is far too large to ever instantiate as a tensor.

`SUNNet` + `sun_value_cx` is the path the codegen actually uses. It contracts the network but never
builds a tensor: it folds structure constants and generator traces *sparsely* over the typed-out
SU(N) tables (`network/sun_net.hpp`) straight to a single complex number. It holds nets a
materialised tensor never could, and the generated kernel carries only the resulting **number**
$C_F$ — never a colour tensor. This tutorial computes $C_F$ that way and checks it against its
known closed form.

## The program

`Tutorials/01-color-contraction/colour_factors.cpp` — folds the network to a number with
`network::sun_value_cx`, and checks it against the closed forms $\mathrm{tr}(T^aT^a) = 4$ and
$f^{abc}f^{abc} = 24$.

```cpp
#include "numtracer/network/sun_net.hpp" // sun_value_cx + SUN::f / SUN::T / SUN::deltaAdj / SUN::deltaFund

#include <cmath>
#include <cstdio>

using namespace numtracer;
using namespace numtracer::network;

// A single unscoped enum keeps every colour label distinct (auto-numbered), so the adjoint and
// fundamental indices never collide (axes contract iff their labels are equal).
enum {
  a, b, c,  // adjoint (gluon) indices
  A, B      // fundamental (quark) indices
};

int main() {
  // SUN::T(g, a, A, B) = (T^a)_{AB} in SU(g). Sharing the gluon index a sums it; the
  // fundamental labels A -> B -> A close the quark line. Closed: tr(T^a T^a) = (N^2-1)/2.
  SUNNet trTT = {SUN::T(3, a, A, B), SUN::T(3, a, B, A)};
  const Cx t = sun_value_cx(trTT); // = 4 for SU(3)
  const double CF = t.re / 3.0;        // C_F = tr / N = 4 / 3

  // f^{abc} f^{abc} = N(N^2-1) = 24 for SU(3). SUN::f(g, a, b, c) = f^{abc}.
  const Cx ff = sun_value_cx({SUN::f(3, a, b, c), SUN::f(3, a, b, c)});

  std::printf("tr(T^a T^a)      = %g   (expect 4)\n", t.re);
  std::printf("C_F = tr / N     = %g   (expect 4/3 = %g)\n", CF, 4.0 / 3.0);
  std::printf("f^{abc} f^{abc}  = %g   (expect 24)\n", ff.re);

  const bool ok = approx(t, Cx{4, 0}) && std::fabs(CF - 4.0 / 3.0) < 1e-12 &&
                  approx(ff, Cx{24, 0});
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
```

```bash
cmake --build build --target colour_factors && ./build/colour_factors
```
```text
tr(T^a T^a)      = 4   (expect 4)
C_F = tr / N     = 1.33333   (expect 4/3 = 1.33333)
f^{abc} f^{abc}  = 24   (expect 24)
ALL TESTS PASSED
```

## Reading it

A colour network is a `SUNNet` — a list of factors, each built by a small helper. Every label
is an integer you choose (here named through one unscoped `enum`); a label that appears twice is
summed.

| factor              | meaning                                                                               |
| ------------------- | ------------------------------------------------------------------------------------- |
| `SUN::T(g, a, i, j)` | generator $(T^a)_{ij}$ in SU($g$): adjoint index `a`, fundamental row `i`, column `j` |
| `SUN::f(g, a, b, c)`   | structure constant $f^{abc}$ in SU($g$)                                               |
| `SUN::deltaAdj(g, a, b)`  | adjoint Kronecker $\delta^{ab}$                                                       |
| `SUN::deltaFund(g, i, j)`   | fundamental Kronecker $\delta_{ij}$                                                   |

In `trTT` the two generators share the **adjoint** label `a` (the gluon, summed) and chain the
**fundamental** labels `A → B → A` (the quark line, closed into a loop). With every index
summed, `sun_value_cx` returns a scalar: $\mathrm{tr}(T^aT^a) = (N^2-1)/2 = 4$, so
$C_F = 4/N = 4/3$ — matching the known closed form without ever materialising a colour tensor.

As foreshadowed above, when NumTracer generates a kernel this sparse fold runs at build time, so the
kernel carries the *number* $C_F$, never a colour tensor. For the tables and the runtime oracle behind
them, see [Sector data](../internals/sectors.md).

A colour factor need not collapse all the way to one number: a $\delta$ tagged with a
per-component dressing array folds instead to a polynomial $\sum_a c_a Z_a$, which is how
NumTracer dresses individual colour/flavour components — see
[per-flavour and per-component dressings](dressed-flavour.md).

Next: the Dirac trace — the part that carries the momenta — [hand-coded](dirac-traces.md).
