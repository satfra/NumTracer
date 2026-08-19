// REGRESSION — a FUNDAMENTAL flavour index with NO flavour generators anywhere (the kaon case).
//
// A flavour label rides along a line as a Kronecker delta. NumTracer has a cheap path (stitch a
// chain, close a loop to Nf) and the SU(N) engine, and it used to pick the engine only when it also
// saw flavour GENERATORS. A kaon-like field (K_i ~ sbar g5 q_i) carries a fundamental flavour index
// but couples through no generator, so it got the cheap path — and in the TADPOLE, where the quartic
// contributes delta_ac delta_bd + delta_ad delta_bc against the loop's own delta, the deltas form a
// WEB rather than a chain. No local rewrite closes a web, so the residue was CForm'd into the kernel
// as `NumTracer_Private_flavDelta(F1, F2)`, an undeclared identifier that both GCC and Clang parse.
//
// THE ORACLE. TFlav[0,f1,f2] = deltaFundFlav[f1,f2]/Sqrt[2 Nf], so the two kernels below are the
// SAME object written two ways: `fund` spells the delta plainly (the path that used to leak),
// `tflav` spells it through the generator route (which always worked, because it trips the generator
// gate — and is the workaround the bug report was using). They must agree pointwise. This is an
// exact self-contained equivalence: no external truth, no dense backend, no FORM oracle.
//
// Note the two kernels are expected to come out not merely equal but BYTE-identical up to their
// namespace, since the fix makes the plain spelling take exactly the generator route. The pointwise
// check is deliberately the weaker statement: it stays meaningful if a future emission change makes
// the two texts diverge while keeping the algebra the same.
#include "Kaon_tadpole_fund_kernel.hh"  // plain deltaFundFlav      (gen/gen_kaon_tadpole.wls)
#include "Kaon_tadpole_tflav_kernel.hh" // Sqrt[2 Nf] TFlav[0,..]   (gen/gen_kaon_tadpole.wls)
#include "shim.hpp"                     // DiFfRG::Fn, ShimRegulator

#include <cmath>
#include <cstdio>
#include <random>

namespace {
// Boson propagator dressing and its regulator-dot insertion. Shapes are arbitrary but smooth,
// positive and distinguishable, so a dropped or swapped factor cannot cancel numerically.
double gk(double x) { return 1.0 / (1.0 + 0.35 * x * x) + 0.60; }
double gkDot(double x) { return 0.45 + 0.20 * std::sin(0.7 * x); }
} // namespace

int main()
{
  using Fund = DiFfRG::Kaon_tadpole_fund_kernel;
  using TFlav = DiFfRG::Kaon_tadpole_tflav_kernel;

  const DiFfRG::Fn GK{&gk}, GKdot{&gkDot};

  std::mt19937_64 rng(20260813);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);

  double maxRel = 0.0, maxAbs = 0.0;
  const int N = 200000;
  for (int i = 0; i < N; ++i) {
    const double l1 = U(rng), c1 = Uc(rng), p = U(rng);
    const double fu = static_cast<double>(Fund::kernel(l1, c1, p, GK, GKdot));
    const double tf = static_cast<double>(TFlav::kernel(l1, c1, p, GK, GKdot));
    maxRel = std::max(maxRel, std::fabs(fu - tf) / (1e-300 + std::fabs(tf)));
    maxAbs = std::max(maxAbs, std::fabs(fu));
  }

  std::printf("[fundamental flavour delta with NO generators: plain spelling vs the TFlav route]\n");
  std::printf("  deltaFundFlav vs Sqrt[2 Nf] TFlav[0,..]   max rel err = %.3e\n", maxRel);
  std::printf("  max |kernel| = %.3e\n", maxAbs);

  // maxAbs guards against the vacuous pass: a kernel that collapsed to 0 would "agree" perfectly
  // and prove nothing, which is exactly the failure mode a delta-web bug could produce.
  const bool ok = maxRel < 1e-10 && maxAbs > 1e-6;
  std::printf("\n%s\n", ok ? "ALL TESTS PASSED" : "TESTS FAILED");
  return ok ? 0 : 1;
}
