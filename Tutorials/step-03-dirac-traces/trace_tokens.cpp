// step-03b — The same Dirac trace, described as a chain of tokens.
//
// step-03a built the gamma matrices by hand. Usually you instead DESCRIBE the closed chain as a list
// of tokens and let the frame trace it — this is how generated code builds every diagram. A token is
// one factor of the trace-ordered chain:
//   nt::slash(k)  — a slashed momentum k/ = gamma . k
//   nt::gamma(mu) — a free gamma^mu whose index mu stays OPEN until a Lorentz network contracts it
// The list is cyclic, so it already closes into a trace.
#include <numtracer.hpp> // the whole NumTracer API

#include <cmath>
#include <cstdio>

namespace nt = numtracer;
using nt::Cx;

int main() {
  // @snip begin: comp
  // Two momenta with fully symbolic components: p = (p0..p3), q = (q0..q3).
  nt::Frame F;
  auto P0 = F.symbol("p0"), P1 = F.symbol("p1"), P2 = F.symbol("p2"), P3 = F.symbol("p3");
  auto Q0 = F.symbol("q0"), Q1 = F.symbol("q1"), Q2 = F.symbol("q2"), Q3 = F.symbol("q3");
  nt::Momentum p = F.momentum(P0, P1, P2, P3);
  nt::Momentum q = F.momentum(Q0, Q1, Q2, Q3);
  auto [mu, nu] = F.indices<2>();
  // @snip end: comp

  // @snip begin: chain-pq
  // (1) tr(p/ q/): a closed chain of two slashed momenta, nothing left to contract.
  nt::Poly trPQ = F.trace({nt::slash(p), nt::slash(q)});
  // @snip end: chain-pq

  // @snip begin: chain-g
  // (2) tr(gamma^mu p/ gamma_mu q/): two free gammas whose indices mu, nu are tied together by the
  // metric delta_{mu nu}, handed to trace() as the chain's Lorentz network.
  nt::Poly trG = F.trace({nt::gamma(mu), nt::slash(p), nt::gamma(nu), nt::slash(q)}, nt::metric(mu, nu));
  // @snip end: chain-g

  // Evaluate both at one point.
  const double x[8] = {1.0, 0.5, -0.3, 0.2, 0.8, -0.4, 1.2, 0.1};
  const nt::Point pt = F.at(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]);
  double pq = 0;
  for (int m = 0; m < 4; ++m) pq += x[m] * x[4 + m];
  const Cx vPQ = F.eval(trPQ, pt);
  const Cx vG = F.eval(trG, pt);

  std::printf("tr(p/ q/)           = %g   (= 4 p.q = %g, %d monomials)\n", vPQ.re, 4.0 * pq, trPQ.size());
  std::printf("tr(g^mu p/ g_mu q/) = %g   (= -2 tr(p/ q/) = %g)\n", vG.re, -8.0 * pq);

  // gamma^mu a/ gamma_mu = -2 a/ in 4 dimensions gives the second relation.
  const bool ok = std::fabs(vPQ.re - 4.0 * pq) < 1e-12 && std::fabs(vG.re + 8.0 * pq) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
