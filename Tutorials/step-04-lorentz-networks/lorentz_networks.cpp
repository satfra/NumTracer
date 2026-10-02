// step-04 — Lorentz networks contract to a scalar polynomial.
//
// A network of metrics, momentum vectors and transverse projectors has no free indices once it is
// contracted: the frame sums every shared Lorentz index away and collapses the network to a
// polynomial in the frame's symbols, with each surviving inverse propagator 1/l^2 carried along as a
// separate "atom". We apply it to the warm-up p.P(l).p and read off p^2 (1 - cos^2 theta) — the
// angular factor of the ghost/gluon loop.
#include <numtracer.hpp> // the whole NumTracer API

#include <cmath>
#include <cstdio>

namespace nt = numtracer;

int main() {
  // @snip begin: frame
  // The loop frame: p along axis 0, l at angle theta to it in the 0-1 plane. The symbols are the
  // magnitudes p, l and cos(theta); sin(theta) is not a separate input — F.angle derives it, and
  // knows cos^2 + sin^2 = 1, so l^2 = l^2 cos^2 + l^2 sin^2 simplifies to the single monomial l^2.
  nt::Frame F;
  auto P = F.symbol("p"), L = F.symbol("l");
  auto [C, S] = F.angle("theta");
  nt::Momentum p = F.momentum(P, 0, 0, 0);
  nt::Momentum l = F.momentum(L * C, L * S, 0, 0);
  auto [mu, nu] = F.indices<2>();
  // @snip end: frame

  // @snip begin: net
  // The network as one product of three factors:
  //   vec(mu, p)       : p_mu
  //   projT(mu, nu, l) : P_T(l)_{mu nu} = delta_{mu nu} - l_mu l_nu / l^2
  //   vec(nu, p)       : p_nu
  // Sharing mu and nu sums both away: a scalar network, p.P(l).p. Its 1/l^2 is registered by the
  // frame (from l's components) and evaluated with it.
  nt::Poly poly = F.contract(nt::vec(mu, p) * nt::projT(mu, nu, l) * nt::vec(nu, p));
  // @snip end: net

  // Evaluate at one point: values of the independent symbols p, l, cos_theta, in declaration order.
  const double pv = 1.3, lv = 0.86, c = 0.58;
  const double val = F.eval(poly, F.at(pv, lv, c)).re;

  std::printf("contracted monomials = %d   (p.P.p = p^2 - (p.l)^2 / l^2)\n", poly.size());
  std::printf("p.P(l).p             = %g   (= p^2 (1 - cos^2) = %g)\n", val, pv * pv * (1 - c * c));
  std::printf("p.P(l).p / p^2       = %g   (= 1 - cos^2 theta = %g)\n", val / (pv * pv), 1 - c * c);

  const bool ok = std::fabs(val - pv * pv * (1 - c * c)) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
