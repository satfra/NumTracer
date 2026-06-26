// Tutorial 4 — Hand-coding a full diagram, end to end.
//
// This is the capstone: the quark self-energy's Dirac-trace numerator,
//   T_num = tr[ p/ gamma^mu q/ gamma^nu ] P_{mu nu}(l),   q = l - p,
// assembled by hand with the numeric engine and validated two ways at one
// frame point:
//   (1) the numeric path  — describe the chain + Lorentz net as tokens,
//   contract numerically; (2) the closed form   — 4 p ( -3 c l1 + p
//   + 2 c^2 p ), the analytic value for q = l - p.
// The two agree, which is what pins the engine's algebra down.
#include "numtracer/numeric/numeric_contract.hpp" // numeric_value, nprojT, dslash, dgamma

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace nm = numtracer::numeric;
namespace net = numtracer::network;
using numtracer::Cx;

// Name the axis labels. A single unscoped enum keeps every label distinct
// (auto-numbered). rho and sigma are the Lorentz dummies summed into the p- and
// q-slashes; mu and nu are the free gamma legs the projector P_{mu nu}(l) ties
// together.
enum {
  rho,
  mu,
  sigma,
  nu // Lorentz indices
};

int main() {
  // One-angle frame: p along axis 0, l at angle theta in the 0-1 plane, q = l -
  // p.
  const double p = 1.3, l0 = 0.5, l1y = 0.7; // |p|, and l = (l0, l1y, 0, 0)
  const double pvec[4] = {p, 0, 0, 0};
  const double lvec[4] = {l0, l1y, 0, 0};
  const double qvec[4] = {l0 - p, l1y, 0, 0}; // q = l - p
  const double l2 = l0 * l0 + l1y * l1y;      // l^2

  // ---- (1) the numeric path
  // ---------------------------------------------------------------- Momenta
  // are plain numbers here, so the components are CONSTANT polynomials (nsym =
  // 0 symbols). comp[vid][m] is component m of momentum vid:  p -> 0, q -> 1, l
  // -> 2.
  const int nsym = 0;
  std::vector<std::array<nm::MPoly, 4>> comp(3);
  auto setv = [&](int vid, const double v[4]) {
    for (int m = 0; m < 4; ++m)
      comp[static_cast<std::size_t>(vid)][static_cast<std::size_t>(m)] =
          nm::MPoly::constant(nsym, Cx{v[m], 0});
  };
  setv(0, pvec);
  setv(1, qvec);
  setv(2, lvec);

  // The closed chain p/ gamma^mu q/ gamma^nu, with the two free gammas on
  // Lorentz legs mu, nu.
  net::DiracNet chain = {net::dslash({{1.0, 0}}), net::dgamma(mu),
                         net::dslash({{1.0, 1}}), net::dgamma(nu)};
  // The Lorentz network is the transverse projector P_{mu nu}(l): legs mu/nu,
  // momentum vid 2, and its 1/l^2 factor is tracked as inverse-atom id 0 (its
  // value supplied in atomDen below).
  nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nprojT(mu, nu, {{1.0, 2}}, 0)}}};
  nm::MPoly l2poly(nsym);
  for (int m = 0; m < 4; ++m)
    l2poly = l2poly + comp[2][static_cast<std::size_t>(m)] *
                          comp[2][static_cast<std::size_t>(m)];
  std::vector<nm::MPoly> atomDen = {
      l2poly}; // atom 0 = l^2 (so its reciprocal is 1/l^2)

  nm::MPoly tr = nm::numeric_value(nsym, chain, lor, comp, atomDen);
  Cx num = nm::eval(tr, /*symbols*/ {}, /*atom values*/ {1.0 / l2});
  const double numeric_val = num.re;

  // ---- (2) the closed form
  // -----------------------------------------------------------------
  const double l1 = std::sqrt(l2), c = l0 / l1;
  const double closed = 4.0 * p * (-3.0 * c * l1 + p + 2.0 * c * c * p);

  std::printf("numeric  T_num = %.6f\n", numeric_val);
  std::printf(
      "closed   T_num = %.6f   (4 p(-3 c l1 + p + 2 c^2 p), q = l - p)\n",
      closed);

  const bool ok = std::fabs(numeric_val - closed) < 1e-10;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
