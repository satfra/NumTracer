// Tutorial 1 — Colour factors fold to a number.
//
// The quark self-energy exchanges one gluon, so its colour structure is T^a_{ij} T^a_{jk} =
// C_F delta_ik with C_F = (N^2-1)/2N = 4/3 for SU(3). A colour network with no free indices is
// just a number: numtracer::network::sun_value_cx contracts it over the typed-out SU(N)
// tables. We compute C_F this way and check it against its known closed form.
#include "numtracer/network/sun_net.hpp" // sun_value_cx + SUN::f / SUN::T / SUN::deltaAdj / SUN::deltaFund

#include <cmath>
#include <cstdio>

using namespace numtracer;
using namespace numtracer::network;

// Name the colour axis labels. A single unscoped enum keeps every label distinct (auto-numbered),
// so adjoint and fundamental indices never collide (axes contract iff their labels are equal).
enum {
  a, b, c,  // adjoint (gluon) indices
  A, B      // fundamental (quark) indices
};

int main() {
  // SUN::T(g, a, A, B) = (T^a)_{AB} in SU(g): adjoint index a, fundamental row A, column B.
  // Sharing the gluon index a sums it; the fundamental labels A -> B -> A close the
  // quark line into a loop. The closed trace is tr(T^a T^a) = (N^2-1)/2 = C_F * N.
  SUNNet trTT = {SUN::T(3, a, A, B), SUN::T(3, a, B, A)};
  const Cx t = sun_value_cx(trTT); // = 4 for SU(3)
  const double CF = t.re / 3.0;        // C_F = tr / N = 4 / 3

  // A second classic, fully closed: f^{abc} f^{abc} = N(N^2-1) = 24 for SU(3).
  // SUN::f(g, a, b, c) = f^{abc}; the two copies share all three adjoint indices a, b, c.
  const Cx ff = sun_value_cx({SUN::f(3, a, b, c), SUN::f(3, a, b, c)});

  std::printf("tr(T^a T^a)      = %g   (expect 4)\n", t.re);
  std::printf("C_F = tr / N     = %g   (expect 4/3 = %g)\n", CF, 4.0 / 3.0);
  std::printf("f^{abc} f^{abc}  = %g   (expect 24)\n", ff.re);

  const bool ok = approx(t, Cx{4, 0}) && std::fabs(CF - 4.0 / 3.0) < 1e-12 &&
                  approx(ff, Cx{24, 0});
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
