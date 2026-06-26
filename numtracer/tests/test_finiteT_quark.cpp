// Finite-temperature quark self-energy — the worked finite-T example, dependency-free.
//
// In the heat-bath rest frame (component 0 = Euclidean time / Matsubara) the quark wave-function
// self-energy carries the closed Dirac trace
//
//     N(p,q,l) = tr[ (γ·p)  γ^μ  (γ·q)  γ^ν ] · G_{μν}(l) ,     q = internal quark momentum,
//
// contracted with the FINITE-T gluon line G_{μν}(l) = ZAE·P_E(l) + ZAM·P_M(l): the electric
// (time-like-transverse) and magnetic (spatial-transverse) projectors the gluon propagator splits
// into once O(4) breaks. The numeric engine folds the Dirac chain (4×4 chiral matrix products)
// against the Lorentz net (nprojE / nprojM, with the spatial 1/|l⃗|² atom). It is checked against the
// closed-form Euclidean trace identity oracle
//
//     N = 4 [ 2 (p · G · q) − (p·q) tr G ] ,   tr[γ^aγ^μγ^bγ^ν] = 4(δ^{aμ}δ^{bν} − δ^{ab}δ^{μν} + δ^{aν}δ^{μb}),
//
// with p·G·q = ZAE (p·P_E·q) + ZAM (p·P_M·q) and tr G = ZAE·tr P_E + ZAM·tr P_M = ZAE + 2 ZAM.
//
// The external quark's temporal component is the lowest fermionic Matsubara frequency, vec[p,0] = πT
// (T random per point); the loop's temporal component l0 (a summed Matsubara frequency) and all
// spatial components are random — a genuinely broken-O(4) Dirac × projector trace. Prints
// ALL TESTS PASSED / exits non-zero on failure.
#include "numtracer/numeric/numeric_contract.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using numtracer::Cx;
namespace nm = numtracer::numeric;
namespace network = numtracer::network;

// 8 symbols: external p = vars 0..3 (vid 0), loop l = vars 4..7 (vid 1). The internal quark momentum
// is q = l − p (a vlc), the gluon projector momentum is l, atoms: 0 = l², 1 = |l⃗|².
constexpr int nsym = 8;
constexpr int pVid = 0, lVid = 1;
constexpr double ZAE = 1.3, ZAM = 0.7; // placeholder transverse/longitudinal gluon dressings

int main()
{
  std::vector<std::array<nm::MPoly, 4>> comp(2);
  for (int mu = 0; mu < 4; ++mu) {
    comp[pVid][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 0 + mu);
    comp[lVid][static_cast<std::size_t>(mu)] = nm::MPoly::var(nsym, 4 + mu);
  }
  // l² (atom 0) and spatial |l⃗|² (atom 1, component 0 = temporal dropped).
  nm::MPoly l2(nsym), ls2(nsym);
  for (int mu = 0; mu < 4; ++mu)
    l2 = l2 + comp[lVid][static_cast<std::size_t>(mu)] * comp[lVid][static_cast<std::size_t>(mu)];
  for (int mu = 1; mu < 4; ++mu)
    ls2 = ls2 + comp[lVid][static_cast<std::size_t>(mu)] * comp[lVid][static_cast<std::size_t>(mu)];
  const std::vector<nm::MPoly> atomDen = {l2, ls2};

  // closed Dirac trace tr[ (γ·p) γ^μ (γ·q) γ^ν ], free legs μ=100, ν=101; q = l − p.
  const network::DiracNet chain = {network::dslash({{1.0, pVid}}), network::dgamma(100),
                                   network::dslash({{1.0, lVid}, {-1.0, pVid}}), network::dgamma(101)};
  // finite-T gluon line G_{μν}(l) = ZAE·P_E(l) + ZAM·P_M(l) on the free legs.
  const nm::NNet gluon = {nm::NTerm{Cx{ZAE, 0}, {nm::nprojE(100, 101, {{1.0, lVid}}, 0, 1)}},
                          nm::NTerm{Cx{ZAM, 0}, {nm::nprojM(100, 101, {{1.0, lVid}}, 1)}}};

  const nm::MPoly N = nm::numeric_value(nsym, chain, gluon, comp, atomDen);

  std::mt19937_64 rng(20260622);
  std::uniform_real_distribution<double> Usp(-2.0, 2.0), UT(0.05, 0.6);
  constexpr double tol = 1e-10, PI = 3.14159265358979323846;
  double worst = 0;
  const int npts = 5000;

  for (int it = 0; it < npts; ++it) {
    std::vector<double> x(nsym);
    const double T = UT(rng);
    x[0] = PI * T;                    // vec[p,0] = πT  (lowest fermionic Matsubara frequency)
    for (int mu = 1; mu < 4; ++mu)
      x[static_cast<std::size_t>(mu)] = Usp(rng); // spatial p⃗
    for (int mu = 0; mu < 4; ++mu)
      x[static_cast<std::size_t>(4 + mu)] = Usp(rng); // loop l (incl. independent temporal l0)

    double pp[4], ll[4], qq[4], l2v = 0, ls2v = 0;
    for (int mu = 0; mu < 4; ++mu) {
      pp[mu] = x[static_cast<std::size_t>(mu)];
      ll[mu] = x[static_cast<std::size_t>(4 + mu)];
      qq[mu] = ll[mu] - pp[mu];
      l2v += ll[mu] * ll[mu];
      if (mu > 0) ls2v += ll[mu] * ll[mu];
    }
    const std::vector<double> atomVal = {1.0 / l2v, 1.0 / ls2v};
    const double num = nm::eval(N, x, atomVal).re;

    // closed-form oracle
    double PE[4][4], PM[4][4];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const double PT = (i == j ? 1.0 : 0.0) - ll[i] * ll[j] / l2v;
        PM[i][j] = (i == j && i > 0 ? 1.0 : 0.0) - ((i > 0 && j > 0) ? ll[i] * ll[j] / ls2v : 0.0);
        PE[i][j] = PT - PM[i][j];
      }
    auto pGq = [&](const double M[4][4]) {
      double s = 0;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          s += pp[i] * M[i][j] * qq[j];
      return s;
    };
    double pdotq = 0;
    for (int mu = 0; mu < 4; ++mu)
      pdotq += pp[mu] * qq[mu];
    const double trG = ZAE * 1.0 + ZAM * 2.0; // tr P_E = 1, tr P_M = 2
    const double oracle = 4.0 * (2.0 * (ZAE * pGq(PE) + ZAM * pGq(PM)) - pdotq * trG);

    worst = std::max(worst, std::abs(num - oracle) / (1e-300 + std::abs(oracle)));
  }

  std::printf("== finite-T quark self-energy (%d random points, vec[p,0]=πT) ==\n", npts);
  std::printf("  numeric engine vs closed-form trace identity:  worst relative = %.3e\n", worst);
  const bool ok = worst < tol;
  std::printf(ok ? "\nALL TESTS PASSED\n" : "\nTESTS FAILED\n");
  return ok ? 0 : 1;
}
