// DEMONSTRATION + validation — dressing the u-quark and d-quark DIFFERENTLY.
//
// The generated kernel (gen/gen_flavour_split.wls) is a scalar two-point self-energy from a quark
// loop, summed over two EXPLICIT flavours u and d, each its own diagram with its own propagator
// dressing (Gu vs Gd) and regulator-dot insertion (GuDot vs GdDot) — NO SU(Nf) flavour group, so
// the two dressings flow through INDEPENDENTLY. The emitted kernel is
//   2 * [ tr0 * Gd(|l1|)^2 Gd(|l1-p|) GdDot(|l1|)   +   tr1 * Gu(|l1|)^2 Gu(|l1-p|) GuDot(|l1|) ]
// i.e. (d-loop with Gd) + (u-loop with Gu). This driver:
//   1. validates the numeric kernel against the Dense DTensor oracle with DISTINCT u/d dressings;
//   2. proves the dressings are INDEPENDENT: kernel = (u-only) + (d-only), isolated by zeroing the
//      other flavour's dressing;
//   3. shows flavour-symmetry breaking is real: with Gu != Gd the result differs substantially
//      from the flavour-symmetric (common-dressing) kernel — it is NOT just Nf x (one dressing).
#include "Flavour_split_num_dense_kernel.hh" // Dense DTensor oracle  (gen/)
#include "Flavour_split_num_kernel.hh"       // numeric kernel        (gen/)
#include "shim.hpp"                           // DiFfRG::Fn, ShimRegulator

#include <cmath>
#include <cstdio>
#include <random>

namespace {
// Distinct, smooth, positive per-flavour dressing functions. u and d differ markedly (a lighter u
// dressing vs a heavier d dressing); the regulator-dot insertions differ too.
double gu(double x) { return 1.0 / (1.0 + 0.2 * x) + 0.85; }
double gd(double x) { return 1.0 / (1.0 + 0.6 * x) + 0.40; }
double guDot(double x) { return 0.50 + 0.10 * std::sin(0.3 * x); }
double gdDot(double x) { return 0.70 + 0.20 * std::cos(0.25 * x); }
double zero(double) { return 0.0; }
} // namespace

int main()
{
  using Dense = DiFfRG::Flavour_split_num_dense_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::Flavour_split_num_kernel<DiFfRG::ShimRegulator>;

  const DiFfRG::Fn Gu{&gu}, Gd{&gd}, GuDot{&guDot}, GdDot{&gdDot}, Z{&zero};

  auto cN = [&](double l1, double c1, double p, const DiFfRG::Fn &u, const DiFfRG::Fn &d,
                const DiFfRG::Fn &ud, const DiFfRG::Fn &dd) {
    return static_cast<double>(Num::kernel(l1, c1, p, u, d, ud, dd));
  };
  auto cD = [&](double l1, double c1, double p, const DiFfRG::Fn &u, const DiFfRG::Fn &d,
                const DiFfRG::Fn &ud, const DiFfRG::Fn &dd) {
    return static_cast<double>(Dense::kernel(l1, c1, p, u, d, ud, dd));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);

  double pwErr = 0.0;   // (1) numeric-vs-dense
  double addErr = 0.0;  // (2) full vs (u-only + d-only)
  double breakRel = 0.0;// (3) full vs flavour-symmetric (common dressing = u's)
  double maxAbs = 0.0;
  const int N = 200000;
  for (int i = 0; i < N; ++i) {
    const double l1 = U(rng), c1 = Uc(rng), p = U(rng);

    const double full = cN(l1, c1, p, Gu, Gd, GuDot, GdDot);
    const double dense = cD(l1, c1, p, Gu, Gd, GuDot, GdDot);
    pwErr = std::max(pwErr, std::fabs(full - dense) / (1e-300 + std::fabs(dense)));
    maxAbs = std::max(maxAbs, std::fabs(full));

    // (2) independence: zero one flavour's dressing to isolate the other loop.
    const double uOnly = cN(l1, c1, p, Gu, Z, GuDot, GdDot); // Gd=0 -> d-term vanishes
    const double dOnly = cN(l1, c1, p, Z, Gd, GuDot, GdDot); // Gu=0 -> u-term vanishes
    addErr = std::max(addErr, std::fabs(full - (uOnly + dOnly)) / (1e-300 + std::fabs(full)));

    // (3) flavour-symmetric reference: both flavours forced to the SAME (u) dressing.
    const double sym = cN(l1, c1, p, Gu, Gu, GuDot, GuDot);
    breakRel = std::max(breakRel, std::fabs(full - sym) / (1e-300 + std::fabs(sym)));
  }

  std::printf("[dressing u-quarks and d-quarks differently — broken flavour symmetry]\n");
  std::printf("  (1) numeric vs dense oracle (distinct Gu,Gd)   max rel err = %.3e\n", pwErr);
  std::printf("  (2) full == u-only + d-only (dressings independent) max rel = %.3e\n", addErr);
  std::printf("  (3) full vs flavour-symmetric (Gu==Gd) max rel diff   = %.3e  (>0 => u/d genuinely differ)\n",
              breakRel);
  std::printf("  max |kernel| = %.3e\n", maxAbs);

  // (1) two evaluation paths of the same kernel: FP round-off only -> documented numeric tol 1e-8.
  // (2) algebraic identity within ONE path -> tight.
  // (3) with Gu != Gd the per-flavour split must produce a clearly different answer.
  const bool ok = pwErr < 1e-8 && addErr < 1e-12 && breakRel > 1e-2 && maxAbs > 1e-6;
  std::printf("\n%s\n", ok ? "ALL TESTS PASSED" : "TESTS FAILED");
  return ok ? 0 : 1;
}
