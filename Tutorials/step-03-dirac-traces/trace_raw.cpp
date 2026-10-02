// step-03a — A Dirac trace, hand-coded as a matrix product.
//
// A trace tr(p/ q/) is literally the trace of a product of 4x4 gamma matrices. NumTracer builds
// each slashed momentum as a 4x4 matrix whose entries are POLYNOMIALS in the momentum components
// (nt::Poly), multiplies the chain, and takes the matrix trace. The identity tr(p/ q/) = 4 p.q falls
// straight out of the product — we check it against that closed form.
#include <numtracer.hpp> // the whole NumTracer API

#include <array>
#include <cmath>
#include <cstdio>

namespace nt = numtracer;
using nt::Cx;

int main() {
  // @snip begin: symbols
  // The 8 components of two momenta p, q are the frame's symbols. A Symbol converts to the
  // polynomial it stands for.
  nt::Frame F;
  const std::array<nt::Symbol, 4> p = {F.symbol("p0"), F.symbol("p1"), F.symbol("p2"), F.symbol("p3")};
  const std::array<nt::Symbol, 4> q = {F.symbol("q0"), F.symbol("q1"), F.symbol("q2"), F.symbol("q3")};
  // @snip end: symbols

  // @snip begin: slash
  // F.slashC(k) = sum_mu k_mu gamma_mu: a slashed momentum as a 4x4 matrix of polynomials.
  nt::Mat4 ps = F.slashC({p[0], p[1], p[2], p[3]});
  nt::Mat4 qs = F.slashC({q[0], q[1], q[2], q[3]});

  // tr(p/ q/) is the trace of the matrix product: a polynomial in the 8 symbols.
  nt::Poly tr = nt::mtrace(nt::matmul(ps, qs));
  // @snip end: slash

  // Evaluate at one point (values in declaration order) and compare with 4 p.q.
  const double x[8] = {1.0, 0.5, -0.3, 0.2, 0.8, -0.4, 1.2, 0.1};
  const Cx got = F.eval(tr, F.at(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]));
  double pq = 0;
  for (int mu = 0; mu < 4; ++mu) pq += x[mu] * x[4 + mu];

  std::printf("tr(p/ q/)  = %g + %gi   (%d monomials in the polynomial)\n", got.re, got.im, tr.size());
  std::printf("4 (p.q)    = %g\n", 4.0 * pq);

  const bool ok = std::fabs(got.re - 4.0 * pq) < 1e-12 && std::fabs(got.im) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
