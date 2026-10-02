// step-18 — Dressed propagator numerators: one trace instead of 2^D diagrams.
//
// A dressed quark propagator numerator is a SUM of Dirac structures whose coefficients are
// *runtime* dressings, e.g.
//
//     S(p) = Mq · 𝟙  +  Z(p) · p̸        (a mass term + a wave-function-dressed slash)
//
// A diagram with D such dressed numerators is, if you distribute the sum, 2^D separate traces —
// one per choice of structure in each numerator. NumTracer does NOT distribute. It keeps each
// numerator EAGER (a "slot" = the list of its structure options) and contracts the chain ONCE,
// collecting the result into a DPoly: a polynomial whose *variables* are the dressing calls and
// whose *coefficients* are the kinematic polynomials (nt::Poly) the engine already computes. The
// dressings ride along as opaque atom-ids and never enter the trace arithmetic — so the Dirac /
// Lorentz work is done a single time no matter how many structures each numerator carries.
//
// This is the Dirac-side analogue of step-17's SU(N) fold (there a group-diagonal δ folds to a
// SUNPoly; here a dressed numerator folds to a DPoly). See dpoly.hpp for the type and
// internals/numeric-engine.md for how DPoly wraps Poly.
//
// The example: a quark bubble  tr( γ^μ · S(p) · γ^ν · S(q) )  closed by the gluon metric δ_{μν},
// with two dressed numerators S(p), S(q). Distributed, that is 2×2 = 4 traces; collected, it is
// ONE DPoly of (at most) 4 dressing monomials. We validate that the collected DPoly, evaluated at
// random kinematics and random dressing values, equals the explicit distributed sum to 1e-10.
#include <numtracer.hpp> // the whole NumTracer API

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace nt = numtracer;
using nt::Cx;

static double cdiff(Cx a, Cx b) { return std::abs(a.re - b.re) + std::abs(a.im - b.im); }

int main() {
  // Two momenta p, q with fully symbolic components.
  nt::Frame F;
  auto P0 = F.symbol("p0"), P1 = F.symbol("p1"), P2 = F.symbol("p2"), P3 = F.symbol("p3");
  auto Q0 = F.symbol("q0"), Q1 = F.symbol("q1"), Q2 = F.symbol("q2"), Q3 = F.symbol("q3");
  nt::Momentum p = F.momentum(P0, P1, P2, P3);
  nt::Momentum q = F.momentum(Q0, Q1, Q2, Q3);
  auto [mu, nu] = F.indices<2>();

  // @snip begin: slots
  // The Lorentz half: the two free gluon legs mu, nu meet through one metric delta_{mu nu}.
  const nt::LorentzNet lor = nt::metric(mu, nu);

  // The dressed numerators. A DSlot is the list of a numerator's structure options; a DSlotOpt is
  //   { coeff , dressing-atom ids , Dirac tokens , extra Lorentz-net factors } .
  // The identity structure 𝟙 is the EMPTY token list — a slot option contributes nothing to the
  // chain but its dressing. A slash is one `slash` token. (The `toks`/`netFacs` pair is general
  // enough to hold an open-leg vertex structure too, e.g. `{gamma(mu)}`.)
  //   S(p): option 0 = identity 𝟙 dressed by atom 0 ("Mq");   option 1 = p̸ dressed by atom 1 ("Z(p)").
  //   S(q): option 0 = identity 𝟙 dressed by atom 0 ("Mq");   option 1 = q̸ dressed by atom 2 ("Z(q)").
  // (Atom 0 is shared: both mass terms are the SAME runtime Mq. Atoms 1,2 are the per-leg Z's.)
  nt::DSlot sP = {nt::DSlotOpt{Cx{1, 0}, {0}, /*𝟙*/ {}, {}}, nt::DSlotOpt{Cx{1, 0}, {1}, {nt::slash(p)}, {}}};
  nt::DSlot sQ = {nt::DSlotOpt{Cx{1, 0}, {0}, /*𝟙*/ {}, {}}, nt::DSlotOpt{Cx{1, 0}, {2}, {nt::slash(q)}, {}}};
  // @snip end: slots

  // @snip begin: collect
  // The dressed Dirac chain, in trace order: a token is either a FIXED factor (dtfix) or a SLOT
  // reference (dtslot i -> the i-th entry of the slot list below).  gamma^mu · S(p) · gamma^nu · S(q).
  std::vector<nt::DChainTok> dchain = {nt::dtfix(nt::gamma(mu)), nt::dtslot(0), nt::dtfix(nt::gamma(nu)),
                                       nt::dtslot(1)};

  // Collect: ONE contraction, no 2^D blowup. dp is the DPoly.
  nt::DPoly dp = F.trace(dchain, {sP, sQ}, lor);
  // @snip end: collect

  // Distributed reference: enumerate the 2×2 structure choices, trace each concrete (undressed)
  // chain the ordinary way, and weight by the product of that choice's dressings.
  auto refTrace = [&](int cp, int cq, const nt::Point &pt) {
    nt::DiracChain c = {nt::gamma(mu)};
    if (cp == 1) c.push_back(nt::slash(p)); // p̸ (else identity: nothing to push)
    c.push_back(nt::gamma(nu));
    if (cq == 1) c.push_back(nt::slash(q)); // q̸
    return F.eval(F.trace(c, lor), pt);
  };

  std::printf("Dressed quark bubble  tr( γ^μ S(p) γ^ν S(q) ) δ_{μν},  S = Mq·𝟙 + Z·slash\n");
  std::printf("  DPoly dressing monomials (one trace collected) : %d   (distributed would be 4 traces)\n",
              dp.size());

  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  double maxerr = 0.0;
  int bad = 0;
  for (int it = 0; it < 5000; ++it) {
    std::vector<double> x(8);
    for (double &v : x) v = U(rng);
    const nt::Point pt = F.at(x);
    // One value per dressing id the slot options reference ({0,1,2}), indexed by id:
    // drVal[0]=Mq (shared by both mass terms), drVal[1]=Z(p), drVal[2]=Z(q).
    std::vector<double> drVal = {U(rng), U(rng), U(rng)};

    // collected: evaluate the single DPoly at the kinematic point with these dressing values
    Cx collected = F.eval(dp, pt, drVal);

    // distributed: Σ over the 4 structure choices of (dressing product) · trace(choice)
    Cx dist{0, 0};
    for (int cp = 0; cp < 2; ++cp)
      for (int cq = 0; cq < 2; ++cq) {
        double w = (cp == 0 ? drVal[0] : drVal[1]) * (cq == 0 ? drVal[0] : drVal[2]);
        Cx tr = refTrace(cp, cq, pt);
        dist = dist + Cx{tr.re * w, tr.im * w};
      }

    double e = cdiff(collected, dist);
    maxerr = std::max(maxerr, e);
    if (e >= 1e-10) ++bad;
  }

  const bool ok = (bad == 0);
  std::printf("  collected == distributed over 5000 random points : worst |Δ| = %.2e  (%d bad)\n", maxerr, bad);
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
