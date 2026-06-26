// Tutorial 2b — The same Dirac trace, via the engine's token API.
//
// Tutorial 2a built the gamma matrices by hand. Usually you instead DESCRIBE
// the closed chain as a list of tokens and let numeric_value contract it — this
// is how the code generator builds every diagram. A token is one factor of the
// trace-ordered chain:
//   dslash({{c, vid}})  : a slashed momentum, components = the linear
//   combination c * (momentum vid) dgamma(id)           : a free gamma^id with
//   an OPEN Lorentz index, to be contracted later
// The chain list is cyclic, so it already closes into a trace.
#include "numtracer/numeric/numeric_contract.hpp" // numeric_value, NNet/NTerm, nmet, dslash, dgamma

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace nm = numtracer::numeric;
namespace net = numtracer::network;
using numtracer::Cx;

// Name the free gamma legs' Lorentz indices instead of writing raw label
// numbers. A single unscoped enum keeps every label distinct (auto-numbered).
enum { mu, nu }; // Lorentz indices

int main() {
  // Momentum components as symbols: p -> 0..3, q -> 4..7. comp[vid][m] is the
  // m-th component of momentum vid, here just the symbol itself.
  const int nsym = 8;
  std::vector<std::array<nm::MPoly, 4>> comp(2);
  for (int m = 0; m < 4; ++m) {
    comp[0][static_cast<std::size_t>(m)] = nm::MPoly::var(nsym, m);
    comp[1][static_cast<std::size_t>(m)] = nm::MPoly::var(nsym, 4 + m);
  }

  // (1) tr(p/ q/): a closed chain of two slashed momenta (no free legs, no
  // Lorentz network).
  net::DiracNet chainPQ = {net::dslash({{1.0, 0}}), net::dslash({{1.0, 1}})};
  nm::MPoly trPQ =
      nm::numeric_value(nsym, chainPQ, /*lorentz*/ {}, comp, /*atomDen*/ {});

  // (2) tr(gamma^mu p/ gamma_mu q/): two FREE gamma legs whose Lorentz indices
  // mu, nu are tied together by a metric. A Lorentz network (NNet) is a sum of
  // terms; here one term, coefficient 1, with a single factor nmet(mu, nu) =
  // delta_{mu nu}.
  net::DiracNet chainG = {net::dgamma(mu), net::dslash({{1.0, 0}}),
                          net::dgamma(nu), net::dslash({{1.0, 1}})};
  nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(mu, nu)}}};
  nm::MPoly trG = nm::numeric_value(nsym, chainG, lor, comp, {});

  // Evaluate both at one point.
  std::vector<double> x = {1.0, 0.5, -0.3, 0.2, 0.8, -0.4, 1.2, 0.1};
  double pq = 0;
  for (int m = 0; m < 4; ++m)
    pq += x[static_cast<std::size_t>(m)] * x[static_cast<std::size_t>(4 + m)];
  Cx vPQ = nm::eval(trPQ, x, {});
  Cx vG = nm::eval(trG, x, {});

  std::printf("tr(p/ q/)           = %g   (= 4 p.q = %g, %zu monomials)\n",
              vPQ.re, 4.0 * pq, trPQ.size());
  std::printf("tr(g^mu p/ g_mu q/) = %g   (= -2 tr(p/ q/) = %g)\n", vG.re,
              -8.0 * pq);

  // The Lorentz identity gamma^mu a/ gamma_mu = -2 a/ in 4D gives the second
  // relation.
  const bool ok = std::fabs(vPQ.re - 4.0 * pq) < 1e-12 &&
                  std::fabs(vG.re + 8.0 * pq) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
