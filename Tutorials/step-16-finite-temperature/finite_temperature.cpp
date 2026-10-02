// step-16 — Finite temperature: the heat bath breaks O(4), and the gluon line splits.
//
// At T > 0 the Euclidean O(4) rotation symmetry breaks to spatial O(3): the heat bath singles out a
// rest frame u^mu = (1,0,0,0), and the ZEROTH (Euclidean-time / Matsubara) component of every
// momentum becomes physically distinct from the spatial three. The metric is unchanged — it is
// still the Euclidean delta — so nothing in the contraction engine needs to change. What changes is
// that a momentum now carries an INDEPENDENT temporal component, and that the transverse projector
// splits into two pieces that must be dressed separately.
//
// This program does two things:
//
//   (A) checks the projector algebra itself — tr P_E = 1, tr P_M = 2, P_E + P_M = P_T — entirely
//       through the engine, so it is a statement about the engine and not about an oracle;
//   (B) runs the running example at finite T: the quark wave-function self-energy
//
//           N(p,q,l) = tr[ p̸ γ^μ q̸ γ^ν ] G_{μν}(l),   q = l - p,
//           G_{μν}(l) = ZAE·P_E(l) + ZAM·P_M(l),
//
//       against the closed-form Euclidean trace identity
//
//           N = 4 [ 2 (p·G·q) - (p·q) tr G ] .
//
// The external quark sits at the lowest fermionic Matsubara frequency, p_0 = πT, and the loop
// carries an independent temporal component l_0 — a genuinely broken-O(4) trace, not a vacuum one
// dressed up. (At a T = 0 vacuum frame every q_0 vanishes and both checks would pass vacuously.)
#include <numtracer.hpp> // the whole NumTracer API

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace nt = numtracer;

constexpr double ZAE = 1.3, ZAM = 0.7; // placeholder electric / magnetic gluon dressings
constexpr double PI = 3.14159265358979323846;

int main() {
  // @snip begin: atoms
  // The frame. Every component of both momenta is an independent symbol — including the TEMPORAL
  // ones, which is the whole point: at finite T, p_0 is not tied to |p⃗|.
  nt::Frame F;
  auto P0 = F.symbol("p0"), P1 = F.symbol("p1"), P2 = F.symbol("p2"), P3 = F.symbol("p3");
  auto L0 = F.symbol("l0"), L1 = F.symbol("l1"), L2 = F.symbol("l2"), L3 = F.symbol("l3");
  nt::Momentum p = F.momentum(P0, P1, P2, P3);
  nt::Momentum l = F.momentum(L0, L1, L2, L3);
  auto [mu, nu] = F.indices<2>();
  // The projectors divide by TWO different denominators: the 4d transverse projector by l^2, the
  // magnetic one by the SPATIAL |l⃗|^2 = l^2 - l_0^2. The frame derives both from l's components
  // (component 0 is the heat-bath direction); mixing them up is the characteristic finite-T bug.
  // @snip end: atoms

  // @snip begin: traces
  // (A) The projector algebra, checked through the engine. Contracting a projector with a metric on
  // the SAME two legs is its trace: delta_{mu nu} P_{mu nu} = tr P. These must come out as the pure
  // numbers 1, 2 and 3 — no kinematics at all — because P_E projects onto one direction (the
  // time-like-transverse one), P_M onto the two spatial-transverse ones, and P_T onto all three.
  const nt::Poly trE = F.contract(nt::projE(mu, nu, l) * nt::metric(mu, nu));
  const nt::Poly trM = F.contract(nt::projM(mu, nu, l) * nt::metric(mu, nu));
  const nt::Poly trT = F.contract(nt::projT(mu, nu, l) * nt::metric(mu, nu));
  // @snip end: traces

  // @snip begin: chain
  // (B) The finite-T self-energy numerator. The Dirac chain is identical to the vacuum one from
  // step-05 — the broken symmetry lives entirely in the FRAME and in the PROJECTORS, not in the
  // gamma algebra. The gluon line is a two-term Lorentz network: one term per projector, each with
  // its own dressing.
  const nt::LorentzNet gluon = ZAE * nt::projE(mu, nu, l) + ZAM * nt::projM(mu, nu, l);
  const nt::Poly N = F.trace({nt::slash(p), nt::gamma(mu), nt::slash(l - p), nt::gamma(nu)}, gluon);
  // @snip end: chain

  bool ok = true;
  std::mt19937_64 rng(20260622);
  std::uniform_real_distribution<double> Usp(-2.0, 2.0), UT(0.05, 0.6);
  double worstTr = 0.0, worstN = 0.0;
  const int npts = 5000;

  for (int it = 0; it < npts; ++it) {
    std::vector<double> x(8);
    const double T = UT(rng);
    x[0] = PI * T;                                    // p_0 = πT, the lowest fermionic Matsubara mode
    for (int c = 1; c < 4; ++c) x[c] = Usp(rng);      // spatial p⃗
    for (int c = 0; c < 4; ++c) x[4 + c] = Usp(rng);  // loop l, independent temporal l_0

    double pp[4], ll[4], qq[4], l2v = 0, ls2v = 0;
    for (int c = 0; c < 4; ++c) {
      pp[c] = x[c];
      ll[c] = x[4 + c];
      qq[c] = ll[c] - pp[c];
      l2v += ll[c] * ll[c];
      if (c > 0) ls2v += ll[c] * ll[c];
    }
    const nt::Point pt = F.at(x);

    // (A) the traces must be constants
    worstTr = std::max({worstTr, std::fabs(F.eval(trE, pt).re - 1.0), std::fabs(F.eval(trM, pt).re - 2.0),
                        std::fabs(F.eval(trT, pt).re - 3.0)});

    // (B) the self-energy against the closed-form trace identity
    double PEm[4][4], PMm[4][4];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const double t = (i == j ? 1.0 : 0.0) - ll[i] * ll[j] / l2v;
        PMm[i][j] = (i == j && i > 0 ? 1.0 : 0.0) - ((i > 0 && j > 0) ? ll[i] * ll[j] / ls2v : 0.0);
        PEm[i][j] = t - PMm[i][j];
      }
    auto pGq = [&](const double M[4][4]) {
      double s = 0;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) s += pp[i] * M[i][j] * qq[j];
      return s;
    };
    double pdotq = 0;
    for (int c = 0; c < 4; ++c) pdotq += pp[c] * qq[c];
    const double trG = ZAE * 1.0 + ZAM * 2.0; // tr P_E = 1, tr P_M = 2
    const double oracle = 4.0 * (2.0 * (ZAE * pGq(PEm) + ZAM * pGq(PMm)) - pdotq * trG);
    const double got = F.eval(N, pt).re;
    worstN = std::max(worstN, std::fabs(got - oracle) / (1e-300 + std::fabs(oracle)));
  }

  std::printf("A. the thermal projector split, contracted through the engine\n");
  std::printf("   tr P_E = 1, tr P_M = 2, tr P_T = 3 : worst |error| over %d points = %.3e\n", npts,
              worstTr);
  std::printf("B. finite-T quark self-energy (p_0 = πT, independent loop l_0)\n");
  std::printf("   engine vs closed-form trace identity : worst relative = %.3e\n", worstN);

  ok = (worstTr < 1e-10) && (worstN < 1e-10);
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
