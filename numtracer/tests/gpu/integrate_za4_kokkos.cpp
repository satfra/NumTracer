// NumTracer — ZA4 flow integrated on GPU via KOKKOS (CUDA backend); the Kokkos twin of
// integrate_za4_gpu.cu. See integrate_za3_kokkos.cpp.
#include "ZA4_gpu_kernel.hh"
#include "integrator_kokkos.hpp"
#include "quadrature.hpp"
#include "ref.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using K = DiFfRG::ZA4_gpu_kernel<DiFfRG::ShimRegulator>;
static constexpr int NV = 4;

struct EvalZA4 {
  DressingSet d;
  double k;
  KOKKOS_INLINE_FUNCTION double operator()(double l1, double c1, double c2, double phi,
                                           double p) const {
    return numtracer::kokkos::real_of(K::kernel(l1, c1, c2, phi, p, k, d.ZA3, d.ZAcbc, d.ZA4,
                                                d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq,
                                                d.Mq));
  }
};

int main(int argc, char **argv) {
  Kokkos::ScopeGuard guard(argc, argv);
  const double k = 1.0, x_extent = 2.0;
  const double pref = std::pow(2.0 * M_PI, -4.0);
  const EvalZA4 ev{DressingSet{}, k};
  DressingSet d;

  std::vector<double> P(64);
  for (size_t i = 0; i < P.size(); ++i)
    P[i] = std::pow(10.0, -2.0 + 3.0 * (double)i / (double)(P.size() - 1));

  std::vector<Axis> ax;
  ax.push_back(gsl_axis(Quad::legendre, 32, 0.0, std::sqrt(x_extent) * k));
  ax.push_back(gsl_axis(Quad::chebyshev2, 8, -1.0, 1.0));
  ax.push_back(gsl_axis(Quad::legendre, 8, -1.0, 1.0));
  ax.push_back(gsl_axis(Quad::legendre, 8, 0.0, 2.0 * M_PI));

  numtracer::kokkos::SegmentedQuadIntegrator<EvalZA4, NV> I;
  for (int a = 0; a < NV; ++a) I.set_axis(a, ax[a].nodes, ax[a].weights);
  auto run = [&]() {
    auto r = I.map(P, pref, ev);
    for (size_t i = 0; i < P.size(); ++i)
      r.integral[i] += std::real(K::constant(P[i], k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc,
                                             d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
    return r;
  };

  bool ok = true;

  // --- validation vs the CPU same-quadrature reference --------------------------------------
  {
    const auto gpu = run();
    const auto cpu = cpu_reference_za4(ax.data(), P, k, pref);
    double rel = 0;
    for (size_t i = 0; i < P.size(); ++i) {
      if (!std::isfinite(gpu.integral[i])) { ok = false; }
      rel = std::max(rel, std::fabs(gpu.integral[i] - cpu[i]) / std::max(std::fabs(cpu[i]), 1e-300));
    }
    ok &= rel < 1e-8;
    std::printf("validation  orders 32x8x8x8, %zu p-points: max rel err = %.3e  %s\n", P.size(), rel,
                rel < 1e-8 ? "PASS" : "FAIL");
  }

  // --- timing at the same production orders ---------------------------------------------------
  {
    for (int w = 0; w < 2; ++w) run();
    std::vector<numtracer::kokkos::MapResult> rs;
    for (int t = 0; t < 5; ++t) rs.push_back(run());
    for (const auto &r : rs)
      for (double v : r.integral)
        if (!std::isfinite(v)) ok = false;
    auto best = std::min_element(rs.begin(), rs.end(), [](const auto &a, const auto &b) {
                  return a.phase1_ms < b.phase1_ms;
                });
    std::printf("production  orders 32x8x8x8, %zu p-points (%lld evals, buffer %.1f MB) [Kokkos]:\n",
                P.size(), best->evals, I.buffer_bytes() / 1048576.0);
    std::printf("            phase1 (integrand) %.3f ms  (%.1f Mevals/s)\n", best->phase1_ms,
                best->evals / best->phase1_ms / 1e3);
    std::printf("            phase2 (reduction) %.3f ms   wall %.2f ms\n", best->phase2_ms,
                best->wall_ms);
  }

  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
