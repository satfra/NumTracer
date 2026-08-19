// NumTracer — the rank-1 electric projector, checked against the dense one through the real
// contraction.
//
// `push_elem_factors` emits ProjE as an outer product of two 4-vectors instead of one dense 4x4
// factor (see its doc comment). The identity is exact, but "the identity is exact" is not the same
// claim as "the engine contracts both forms to the same polynomial": the two forms take different
// elimination orders, carry different atoms in the numerator, and meet `reduce_units` and the
// monomial-cancellation logic differently. This test makes the second claim directly — contract the
// SAME network both ways, in one process, and compare the resulting MPoly by evaluation.
//
// Comparing by evaluation rather than term-by-term is deliberate: the two forms legitimately produce
// different (algebraically equal) monomial sets, so a structural diff would report differences that
// are not errors. Grading is against the scale of the terms involved, never pointwise-relative — a
// structurally-zero contraction leaves fp residue in one form and an exact 0 in the other, and a
// relative test would call that a 100% error.
#include "numtracer/numeric/env.hpp"
#include "numtracer/numeric/numeric_contract.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace numtracer;
using namespace numtracer::numeric;

namespace
{
  int fail = 0;

  void ok(const char *what, bool cond)
  {
    std::printf("%-64s %s\n", what, cond ? "ok" : "FAIL");
    if (!cond) ++fail;
  }

  // TWO independent momenta: k (vid 0), which the projectors are built from, and p (vid 1), which
  // probes them. Probing with k itself is the trap this test fell into first: P_E is transverse to
  // k, so every k·P_E·k contraction is STRUCTURALLY ZERO, and comparing two zeros — one an exact 0,
  // the other fp residue — measures nothing while looking like a catastrophic failure.
  constexpr int kNsym = 8; // k0..k3, p0..p3

  std::vector<std::array<MPoly, 4>> makeComp()
  {
    LorentzEnv env(kNsym, {});
    std::array<MPoly, 4> k{env.zero(), env.zero(), env.zero(), env.zero()};
    std::array<MPoly, 4> p{env.zero(), env.zero(), env.zero(), env.zero()};
    for (int i = 0; i < 4; ++i) {
      std::vector<int> ek(kNsym, 0), ep(kNsym, 0);
      ek[static_cast<std::size_t>(i)] = 1;
      ep[static_cast<std::size_t>(i) + 4] = 1;
      k[static_cast<std::size_t>(i)] = env.mono(ek, Cx{1., 0});
      p[static_cast<std::size_t>(i)] = env.mono(ep, Cx{1., 0});
    }
    return {k, p};
  }

  /// atomDen[0] = k², atomDen[1] = |k⃗|² — what collect_atom_denoms fills for a ProjE with
  /// inv=0, invS=1. Both are built from k (vid 0) only; p carries no atom.
  std::vector<MPoly> makeAtomDen()
  {
    LorentzEnv env(kNsym, {});
    MPoly k2 = env.zero(), ks2 = env.zero();
    for (int i = 0; i < 4; ++i) {
      std::vector<int> e(kNsym, 0);
      e[static_cast<std::size_t>(i)] = 2;
      const MPoly sq = env.mono(e, Cx{1., 0});
      k2 = k2 + sq;
      if (i > 0) ks2 = ks2 + sq;
    }
    return {k2, ks2};
  }

  /// Evaluate an MPoly at a numeric point. atomVal[a] = 1/atomDen[a].
  double evalAt(const MPoly &poly, const std::array<double, 8> &x8)
  {
    const double k2 = x8[0] * x8[0] + x8[1] * x8[1] + x8[2] * x8[2] + x8[3] * x8[3];
    const double ks2 = x8[1] * x8[1] + x8[2] * x8[2] + x8[3] * x8[3];
    const std::vector<double> x(x8.begin(), x8.end());
    const std::vector<double> atomVal = {1.0 / k2, 1.0 / ks2};
    return eval(poly, x, atomVal).re;
  }

  /// Contract one Lorentz network, forcing either the dense or the rank-1 ProjE factorisation.
  /// Calls the two builders
  /// directly — this test needs both forms in the SAME process.
  MPoly contractBoth(const std::vector<NElem> &elems, bool rank1, const std::vector<MPoly> &aden,
                     const std::vector<std::array<MPoly, 4>> &comp)
  {
    std::vector<ndetail::Factor> facs;
    for (const NElem &el : elems) {
      if (rank1)
        ndetail::push_elem_factors(facs, kNsym, el, comp, aden);
      else
        facs.push_back(ndetail::elem_factor(kNsym, el, comp));
    }
    return ndetail::contract_factors(kNsym, std::move(facs), aden, {});
  }

  struct Gap {
    double rel;    ///< worst |dense - rank1| / (largest value seen in the sweep)
    double maxAbs; ///< largest |value| either form produced — near 0 means structurally zero
    std::size_t denseTerms, rank1Terms;
  };

  /// Contract `elems` both ways over many random momenta and compare.
  ///
  /// A network with no free momentum (a closed trace of orthogonal projectors) is ANALYTICALLY
  /// ZERO, and then `rel` is a ratio of two roundings and means nothing — one form leaves fp
  /// residue, the other cancels exactly. `maxAbs` is what to assert on in that case, which is why
  /// it is returned rather than folded into a single verdict here.
  Gap worstGap(const char *what, const std::vector<NElem> &elems)
  {
    const auto comp = makeComp();
    const auto aden = makeAtomDen();
    const MPoly a = contractBoth(elems, false, aden, comp);
    const MPoly b = contractBoth(elems, true, aden, comp);

    std::mt19937 rng(20260811);
    std::uniform_real_distribution<double> U(-2.0, 2.0);
    double scale = 0., worst = 0.;
    std::vector<std::array<double, 8>> pts;
    for (int t = 0; t < 400; ++t) {
      std::array<double, 8> k{U(rng), U(rng), U(rng), U(rng), U(rng), U(rng), U(rng), U(rng)};
      if (std::fabs(k[1]) + std::fabs(k[2]) + std::fabs(k[3]) < 1e-3) continue; // |k⃗|² ~ 0
      pts.push_back(k);
      scale = std::fmax(scale, std::fabs(evalAt(a, k)));
      scale = std::fmax(scale, std::fabs(evalAt(b, k)));
    }
    const double maxAbs = scale;
    if (scale < 1e-300) scale = 1.0; // both identically zero: avoid 0/0
    for (const auto &k : pts)
      worst = std::fmax(worst, std::fabs(evalAt(a, k) - evalAt(b, k)) / scale);
    std::printf("   %-38s dense=%-5zu rank1=%-5zu terms   rel %.2e   max|v| %.2e\n", what,
                a.terms.size(), b.terms.size(), worst, maxAbs);
    return {worst, maxAbs, a.terms.size(), b.terms.size()};
  }

} // namespace

int main()
{
  // Lorentz index ids: 0..3. inv=0 (k²), invS=1 (|k⃗|²), momentum vid 0.
  const auto E = [](int a, int b) { return nprojE(a, b, {{1.0, 0}}, 0, 1); };
  const auto Tp = [](int a, int b) { return nprojT(a, b, {{1.0, 0}}, 0); };
  const auto M = [](int a, int b) { return nprojM(a, b, {{1.0, 0}}, 1); };
  const auto V = [](int a) { return nvec(a, {{1.0, 1}}); }; // p, NOT k — see makeComp
  const auto G = [](int a, int b) { return nmet(a, b); };

  std::printf("rank-1 vs dense electric projector, through contract_factors:\n");

  // The shapes that matter: a closed trace, contraction against vectors and metrics, and the
  // multi-projector chains a real 4-point finite-T diagram builds.
  const Gap g1 = worstGap("tr P_E", {E(0, 0)});
  const Gap g2 = worstGap("p . P_E . p", {V(0), E(0, 1), V(1)});
  const Gap g3 = worstGap("P_E . P_T (open)", {E(0, 1), Tp(1, 2), V(0), V(2)});
  const Gap g4 = worstGap("tr(P_E P_T)", {E(0, 1), Tp(1, 0)});
  const Gap g5 = worstGap("tr(P_E P_M) == 0 (orthogonal)", {E(0, 1), M(1, 0)});
  const Gap g6 = worstGap("delta . P_E . delta", {G(0, 1), E(1, 2), G(2, 0)});
  const Gap g7 = worstGap("two P_E in one chain", {E(0, 1), Tp(1, 2), E(2, 3), V(0), V(3)});
  const Gap g8 = worstGap("P_E chain x 3 (box-like)",
                             {E(0, 1), Tp(1, 2), E(2, 3), Tp(3, 4), E(4, 5), V(0), V(5)});

  std::printf("\n");
  const double tol = 1e-12;
  ok("tr P_E agrees", g1.rel < tol);
  ok("p.P_E.p agrees", g2.rel < tol);
  ok("P_E.P_T with open legs agrees", g3.rel < tol);
  ok("tr(P_E P_T) agrees", g4.rel < tol);
  // Orthogonality: P_E . P_M is the zero tensor, so BOTH forms must vanish. Asserted absolutely —
  // comparing two roundings of zero to each other is the trap this test is built to avoid.
  ok("tr(P_E P_M) vanishes in both forms", g5.maxAbs < 1e-10);
  ok("metric-sandwiched P_E agrees", g6.rel < tol);
  ok("two electric projectors in a chain agree", g7.rel < tol);
  ok("three electric projectors in a chain agree", g8.rel < tol);
  // The POINT of the change, not just its safety: cutting the network at each electric projector
  // has to shrink the contraction, and by a lot on the multi-projector chains a real finite-T
  // 4-point diagram is made of. If this ever regresses, the rank-1 path is costing correctness
  // review for nothing and should be reverted, not kept.
  ok("rank-1 shrinks the 2-projector chain >=2x", g7.rank1Terms * 2 <= g7.denseTerms);
  ok("rank-1 shrinks the 3-projector chain >=4x", g8.rank1Terms * 4 <= g8.denseTerms);

  std::printf("\n%s (%d failure%s)\n", fail ? "TESTS FAILED" : "ALL TESTS PASSED", fail,
              fail == 1 ? "" : "s");
  return fail ? 1 : 0;
}
