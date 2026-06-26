// DEMONSTRATION + validation — dressing u and d quarks differently WITHIN the SU(2) flavour group.
//
// The in-group kernel (gen/gen_flavour_ingroup.wls) is a SINGLE quark loop carrying an SU(2) flavour
// index whose propagators carry a flavour-DIAGONAL dressing diag(G[0],G[1]); the flavour trace folds
// (via sun_value_dressed) to Σ_{f∈{u,d}} of the per-flavour dressing product, with the Dirac/colour
// trace computed ONCE. It must reproduce the separate-per-flavour-diagram kernel flow_flavour_split
// (gen/gen_flavour_split.wls) EXACTLY: feeding G = {gu, gd}, Gdot = {guDot, gdDot} reproduces the
// (Gu/Gd, GuDot/GdDot) split kernel. This is an exact equivalence oracle (no dense backend needed):
// the two completely different codegen paths must agree to round-off.
#include "Flavour_ingroup_num_kernel.hh" // in-group (SU(2)-diagonal dressing)   (gen/)
#include "Flavour_split_num_kernel.hh"   // separate per-flavour diagrams         (gen/)
#include "shim.hpp"                        // DiFfRG::Fn, ShimRegulator

#include <array>
#include <cmath>
#include <cstdio>
#include <random>

namespace {
// The SAME per-flavour dressing functions used by compare_flavour_split.cpp.
double gu(double x) { return 1.0 / (1.0 + 0.2 * x) + 0.85; }
double gd(double x) { return 1.0 / (1.0 + 0.6 * x) + 0.40; }
double guDot(double x) { return 0.50 + 0.10 * std::sin(0.3 * x); }
double gdDot(double x) { return 0.70 + 0.20 * std::cos(0.25 * x); }
} // namespace

int main()
{
  using InGroup = DiFfRG::Flavour_ingroup_num_kernel<DiFfRG::ShimRegulator>;
  using Split = DiFfRG::Flavour_split_num_kernel<DiFfRG::ShimRegulator>;

  const DiFfRG::Fn Gu{&gu}, Gd{&gd}, GuDot{&guDot}, GdDot{&gdDot};
  // component 0 = u, component 1 = d — matches ntSUNDiagFund's diag(name[0], name[1]).
  const std::array<DiFfRG::Fn, 2> G{Gu, Gd}, Gdot{GuDot, GdDot};

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);

  double maxRel = 0.0, maxAbs = 0.0;
  const int N = 200000;
  for (int i = 0; i < N; ++i) {
    const double l1 = U(rng), c1 = Uc(rng), p = U(rng);
    const double ig = static_cast<double>(InGroup::kernel(l1, c1, p, G, Gdot));
    const double sp = static_cast<double>(Split::kernel(l1, c1, p, Gu, Gd, GuDot, GdDot));
    maxRel = std::max(maxRel, std::fabs(ig - sp) / (1e-300 + std::fabs(sp)));
    maxAbs = std::max(maxAbs, std::fabs(ig));
  }

  std::printf("[u/d dressed differently WITHIN the SU(2) flavour group vs the per-flavour split]\n");
  std::printf("  in-group (diag dressing) vs separate-diagram split   max rel err = %.3e\n", maxRel);
  std::printf("  max |kernel| = %.3e\n", maxAbs);

  // Two unrelated codegen paths (a folded SU(2) trace with a runtime colour-sum token vs two
  // separate diagrams) -> algebraically identical; only FP evaluation-order round-off differs
  // (legal FMA/reassociation regrouping of the 6·Σ vs 2·(Nc·…) forms). Well inside the documented
  // numeric tolerance 1e-8.
  const bool ok = maxRel < 1e-10 && maxAbs > 1e-6;
  std::printf("\n%s\n", ok ? "ALL TESTS PASSED" : "TESTS FAILED");
  return ok ? 0 : 1;
}
