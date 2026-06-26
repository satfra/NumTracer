// Tutorial 3 — Lorentz networks contract to a scalar polynomial.
//
// A network of metrics, momentum vectors and transverse projectors has no "components" once you
// contract it: the numeric engine sums the shared Lorentz indices away over the loop frame and
// collapses the network to a polynomial in the frame's momentum components, carrying any surviving
// inverse propagators 1/l^2 as separate "atoms". numtracer::numeric::numeric_value does that
// contraction. We apply it to the warm-up p.P(l).p and read off p^2 (1 - cos^2 theta) — the angular
// factor of the ghost/gluon loop.
#include "numtracer/numeric/numeric_contract.hpp" // nvec / nprojT (NNet) + numeric_value -> MPoly

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace nm = numtracer::numeric;
namespace network = numtracer::network;
using numtracer::Cx;

// Name the Lorentz axis labels with a single unscoped enum (every label auto-numbered, so all
// distinct); the momentum ids (p = 0, l = 1, indexing the component table) stay plain integers.
enum { mu, nu }; // Lorentz indices

int main() {
  // The loop frame: each momentum is a vid into a component table comp[vid][0..3], whose entries are
  // MPolys. We pick the concrete one-angle frame up front — p along axis 0, l in the 0-1 plane — and
  // keep only its NON-ZERO components symbolic (one scalar variable each), so numeric_value contracts
  // the network into the polynomial in exactly those frame scalars:
  //   p_0 = var 0,   l_0 = var 1,   l_1 = var 2.
  const int nsym = 3;
  std::vector<std::array<nm::MPoly, 4>> comp(2, {nm::MPoly(nsym), nm::MPoly(nsym), nm::MPoly(nsym), nm::MPoly(nsym)});
  comp[0][0] = nm::MPoly::var(nsym, 0); // p = (p_0, 0, 0, 0)
  comp[1][0] = nm::MPoly::var(nsym, 1); // l = (l_0, l_1, 0, 0)
  comp[1][1] = nm::MPoly::var(nsym, 2);

  // The transverse projector P(l)_{mu nu} = delta_{mu nu} - l_mu l_nu / l^2 carries its 1/l^2 as
  // "atom" id 0; numeric_value needs that atom's denominator l^2 = sum_i comp[1][i]^2.
  nm::MPoly l2(nsym);
  for (int i = 0; i < 4; ++i) l2 = l2 + comp[1][i] * comp[1][i];
  std::vector<nm::MPoly> atomDen = {l2};

  // The Lorentz network as one product term (coefficient 1) of three factors:
  //   nvec(Lbl, vlc)           : a 4-vector on Lorentz index Lbl, momentum = sum coeff*comp(vid).
  //   nprojT(Mu, Nu, vlc, atom) : P_{Mu Nu}(l) = delta_{Mu Nu} - l_Mu l_Nu / l^2, l = sum coeff*comp(vid).
  // Sharing index mu between the first p and the projector, and index nu between the projector and
  // the second p, sums both Lorentz indices away -> a scalar network: p.P(l).p.
  nm::NNet net = {nm::NTerm{Cx{1, 0}, {nm::nvec(mu, {{1.0, 0}}),
                                       nm::nprojT(mu, nu, {{1.0, 1}}, 0),
                                       nm::nvec(nu, {{1.0, 0}})}}};

  // Contract: empty Dirac chain (this is a pure-Lorentz network). The result is one MPoly =
  // p.P(l).p = sp(p,p) - sp(p,l)^2 / l^2 = p_0^2 - p_0^2 l_0^2 / l^2 — two monomials in this frame
  // (the second carries the 1/l^2 atom).
  nm::MPoly poly = nm::numeric_value(nsym, network::DiracNet{}, net, comp, atomDen);

  // Evaluate the contracted polynomial at concrete values: p along axis 0, l at angle theta.
  const double Pm = 1.3, l0 = 0.5, l1 = 0.7;
  const std::vector<double> x = {Pm, l0, l1}; // var 0 = p_0, var 1 = l_0, var 2 = l_1
  const double l2v = l0 * l0 + l1 * l1;
  const std::vector<double> atomVal = {1.0 / l2v}; // value of atom 0 = 1/l^2 at this frame
  const double val = nm::eval(poly, x, atomVal).re;

  const double cth = l0 / std::sqrt(l2v); // cos theta
  std::printf("contracted monomials = %zu   (p.P.p = sp(p,p) - sp(p,l)^2 / l^2)\n", poly.size());
  std::printf("p.P(l).p             = %g   (= p^2 (1 - cos^2) = %g)\n", val, Pm * Pm * (1 - cth * cth));
  std::printf("p.P(l).p / p^2       = %g   (= 1 - cos^2 theta = %g)\n", val / (Pm * Pm), 1 - cth * cth);

  const bool ok = std::fabs(val - Pm * Pm * (1 - cth * cth)) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
