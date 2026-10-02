// step-02 — SU(N) factors fold to a number.
//
// The quark self-energy exchanges one gluon, so its colour structure is T^a_{ij} T^a_{jk} =
// C_F delta_ik with C_F = (N^2-1)/2N = 4/3 for SU(3). A closed SU(N) network is just a number, and
// NumTracer folds it exactly. We compute C_F this way and check it, and f^{abc} f^{abc}, against
// their closed forms.
#include <numtracer.hpp> // the whole NumTracer API

#include <cmath>
#include <cstdio>

namespace nt = numtracer;
using nt::Cx;

int main() {
  // @snip begin: labels
  // An SU(N) group object: here colour SU(3). It hands out index labels, and they come in two
  // types — adjoint (a, b, c: extent N^2-1) and fundamental (i, j: extent N) — so an adjoint label
  // cannot end up in a fundamental slot. A flavour SU(2) would be a second object, nt::SUN su2(2).
  nt::SUN su3(3);
  auto [a, b, c] = su3.adjoint<3>();
  auto [i, j] = su3.fundamental<2>();
  // @snip end: labels

  // @snip begin: cf
  // su3.T(a, i, j) = (T^a)_{ij}: adjoint index a, fundamental row i, column j. Sharing a sums the
  // gluon index; i -> j -> i closes the quark line into a loop. The closed trace is
  // tr(T^a T^a) = (N^2-1)/2 = C_F * N.
  const Cx t = su3.value(su3.T(a, i, j) * su3.T(a, j, i)); // = 4 for SU(3)
  const double CF = t.re / 3.0;                            // C_F = tr / N = 4/3
  // @snip end: cf

  // @snip begin: ff
  // A second classic, fully closed: f^{abc} f^{abc} = N(N^2-1) = 24 for SU(3).
  const Cx ff = su3.value(su3.f(a, b, c) * su3.f(a, b, c));
  // @snip end: ff

  std::printf("tr(T^a T^a)      = %g   (expect 4)\n", t.re);
  std::printf("C_F = tr / N     = %g   (expect 4/3 = %g)\n", CF, 4.0 / 3.0);
  std::printf("f^{abc} f^{abc}  = %g   (expect 24)\n", ff.re);

  const bool ok = nt::approx(t, Cx{4, 0}) && std::fabs(CF - 4.0 / 3.0) < 1e-12 && nt::approx(ff, Cx{24, 0});
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
