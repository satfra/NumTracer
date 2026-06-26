// NumTracer — ZA3 flow integrated on GPU: the full 4D loop integral (radial l1 + two angles;
// the trivial phi gives the 2π in the prefactor) for a log-spaced grid of external momenta p,
// in ONE phase-1 launch + one segmented reduction (numtracer/cuda/integrator.cuh, modeled on
// DiFfRG's Integrator_p2_4D_2ang::map). Quadrature: GSL Gauss-Legendre (l1, cos2) +
// Gauss-Chebyshev-2 (cos1, jacobian baked into the weights).
//
//   validation: GPU vs an independently compiled CPU same-quadrature reference (the HOST
//     kernel variant, ref_za3.cpp) at the production orders, all grid points. Tolerance 1e-8:
//     a mechanics bug (index decode, weights, axis order) is O(1), while legitimate deviations
//     (GPU FMA contraction, different summation order) sit around 1e-13.
//   timing: same orders, 2 warm-ups + 5 timed map() calls, per-phase cudaEvent times.
#include "ZA3_gpu_kernel.hh"
#include "numtracer/cuda/integrator.cuh"
#include "quadrature.hpp"
#include "ref.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using K = DiFfRG::ZA3_gpu_kernel<DiFfRG::ShimRegulator>;
static constexpr int NV = 3;

int main() {
  const double k = 1.0, x_extent = 2.0;                      // radial l1 ∈ [0, √x_extent·k]
  const double pref = 2.0 * M_PI / std::pow(2.0 * M_PI, 4);  // phi integral / fourier factor
  DressingSet d;

  // 64 log-spaced external momenta in [1e-2, 1e1]
  std::vector<double> P(64);
  for (size_t i = 0; i < P.size(); ++i)
    P[i] = std::pow(10.0, -2.0 + 3.0 * (double)i / (double)(P.size() - 1));

  auto make_axes = [&](int n0, int n1, int n2) {
    std::vector<Axis> ax;
    ax.push_back(gsl_axis(Quad::legendre, n0, 0.0, std::sqrt(x_extent) * k));
    ax.push_back(gsl_axis(Quad::chebyshev2, n1, -1.0, 1.0));
    ax.push_back(gsl_axis(Quad::legendre, n2, -1.0, 1.0));
    return ax;
  };
  numtracer::cuda::SegmentedQuadIntegrator<K, NV> I;
  auto run = [&](const std::vector<Axis> &ax) {
    for (int a = 0; a < NV; ++a) I.set_axis(a, ax[a].nodes, ax[a].weights);
    return I.map(P, k, pref, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq,
                 d.Mq);
  };

  bool ok = true;

  // production quadrature orders: radial 32, each angle 8 (the typical production sizes used
  // with DiFfRG; radial needs the resolution — regulator structure lives in l1 — while the
  // angular dependence is smooth). Validation runs at the SAME orders, so it checks exactly
  // the production configuration.
  const auto ax = make_axes(32, 8, 8);

  // --- validation vs the CPU same-quadrature reference --------------------------------------
  {
    const auto gpu = run(ax);
    const auto t0 = std::chrono::steady_clock::now();
    const auto cpu = cpu_reference_za3(ax.data(), P, k, pref);
    const double cpu_ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    double rel = 0;
    for (size_t i = 0; i < P.size(); ++i) {
      if (!std::isfinite(gpu.integral[i])) { ok = false; }
      rel = std::max(rel, std::fabs(gpu.integral[i] - cpu[i]) / std::max(std::fabs(cpu[i]), 1e-300));
    }
    ok &= rel < 1e-8;
    std::printf("validation  orders 32x8x8, %zu p-points: max rel err = %.3e  %s\n", P.size(), rel,
                rel < 1e-8 ? "PASS" : "FAIL");
    std::printf("            CPU reference %.1f ms vs GPU map %.2f ms (same %lld evals)\n", cpu_ms,
                gpu.wall_ms, gpu.evals);
    std::printf("            I(p) samples: p=%.3g -> %.6e   p=%.3g -> %.6e   p=%.3g -> %.6e\n", P[0],
                gpu.integral[0], P[P.size() / 2], gpu.integral[P.size() / 2], P.back(),
                gpu.integral.back());
  }

  // --- timing at the same production orders ---------------------------------------------------
  {
    for (int w = 0; w < 2; ++w) run(ax); // warm-up (allocations, module load)
    std::vector<numtracer::cuda::MapResult> rs;
    for (int t = 0; t < 5; ++t) rs.push_back(run(ax));
    for (const auto &r : rs)
      for (double v : r.integral)
        if (!std::isfinite(v)) ok = false;
    auto best = std::min_element(rs.begin(), rs.end(), [](const auto &a, const auto &b) {
                  return a.phase1_ms < b.phase1_ms;
                });
    std::printf("production  orders 32x8x8, %zu p-points (%lld evals, buffer %.1f MB):\n", P.size(),
                best->evals, I.buffer_bytes() / 1048576.0);
    std::printf("            phase1 (integrand) %.3f ms  (%.1f Mevals/s)\n", best->phase1_ms,
                best->evals / best->phase1_ms / 1e3);
    std::printf("            phase2 (reduction) %.3f ms   wall %.2f ms\n", best->phase2_ms,
                best->wall_ms);
  }

  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
