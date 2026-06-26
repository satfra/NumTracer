// Tutorial 2a — A Dirac trace, hand-coded as a matrix product.
//
// A trace tr(p-slash q-slash) is literally the trace of a product of 4x4 gamma matrices. The
// numeric engine builds each slashed momentum as a 4x4 matrix whose entries are POLYNOMIALS in
// the momentum components (numtracer::numeric::MPoly), multiplies the chain, and takes the
// matrix trace. The classic identity tr(p/ q/) = 4 p.q falls straight out of the product — we
// check it against that closed form.
#include "numtracer/numeric/spinor_mat.hpp" // Mat4, slashC, matmul, mtrace (pulls in MPoly, eval)

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace nm = numtracer::numeric;
using numtracer::Cx;

int main() {
  // We treat the 8 components of two momenta p, q as symbols 0..7:
  //   MPoly::var(nsym, i) is "the i-th symbol". p_mu = symbol mu, q_mu = symbol 4+mu.
  const int nsym = 8;
  std::array<nm::MPoly, 4> p, q;
  for (int mu = 0; mu < 4; ++mu) {
    p[static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, mu);
    q[static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 + mu);
  }

  // slashC(nsym, comp) = sum_mu comp[mu] * gamma^mu : the slashed momentum as a 4x4 of MPoly.
  nm::Mat4 ps = nm::slashC(nsym, p);
  nm::Mat4 qs = nm::slashC(nsym, q);

  // tr(p/ q/) is the trace of the matrix product. The result is a polynomial in the 8 symbols.
  nm::MPoly tr = nm::mtrace(nm::matmul(ps, qs));

  // Evaluate at one concrete point and compare to the closed form 4 p.q.
  std::vector<double> x = {1.0, 0.5, -0.3, 0.2, 0.8, -0.4, 1.2, 0.1};
  Cx got = nm::eval(tr, x, /*no inverse atoms*/ {});

  double pq = 0;
  for (int mu = 0; mu < 4; ++mu) pq += x[static_cast<std::size_t>(mu)] * x[static_cast<std::size_t>(4 + mu)];

  std::printf("tr(p/ q/)  = %g + %gi   (%zu monomials in the polynomial)\n", got.re, got.im, tr.size());
  std::printf("4 (p.q)    = %g\n", 4.0 * pq);

  const bool ok = std::fabs(got.re - 4.0 * pq) < 1e-12 && std::fabs(got.im) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
