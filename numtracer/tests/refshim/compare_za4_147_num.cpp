// Cross-check the dense 1/4/7 numeric kernel for ZA4 (four-gluon box, quark loop carrying the full
// AqbqDirect147 vertex) against the FORM-validated struct-1 numeric kernel. No FORM oracle exists at
// this size (rank-4, spFrame3); validation is by STRUCT-REDUCTION: with ZAqbq4 = ZAqbq7 = 0 the 1/4/7
// vertex collapses to structure 1, so the 1/4/7 kernel must reduce EXACTLY (round-off) to
// ZA4_num_kernel (validated vs FORM by compare_za4_num). With the real ZAqbq4/ZAqbq7 it must DIFFER
// (structures 4/7 genuinely contribute). Same precedent as compare_za3_147_num.
#include "ZA4_147_num_kernel.hh" // dense 1/4/7 numeric kernel (gen/), σ^{μν} folded
#include "ZA4_num_kernel.hh"     // struct-1 numeric kernel, validated vs FORM (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using K147 = DiFfRG::ZA4_147_num_kernel<DiFfRG::ShimRegulator>;
  using K1 = DiFfRG::ZA4_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  DiFfRG::Fn zero{+[](double) { return 0.0; }};

  // 1/4/7 kernel with structures 4,7 turned off:
  auto c147_0 = [&](double l, double c1, double c2, double ph, double p, double k) {
    return std::real(K147::kernel(l, c1, c2, ph, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, zero, zero,
                                  d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  // 1/4/7 kernel with the real structure-4,7 dressings:
  auto c147 = [&](double l, double c1, double c2, double ph, double p, double k) {
    return std::real(K147::kernel(l, c1, c2, ph, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4,
                                  d.ZAqbq7, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  // struct-1 kernel (the FORM-validated reference):
  auto c1 = [&](double l, double cc1, double cc2, double ph, double p, double k) {
    return std::real(K1::kernel(l, cc1, cc2, ph, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1,
                                d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(20240619);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999), Uphi(0.0, 6.283185307179586);
  const int Np = 200000;
  double reduceErr = 0, full478Diff = 0;
  for (int i = 0; i < Np; ++i) {
    double l = U(rng), c1v = Uc(rng), c2v = Uc(rng), ph = Uphi(rng), p = U(rng), k = U(rng);
    double ref = c1(l, c1v, c2v, ph, p, k);
    double red = c147_0(l, c1v, c2v, ph, p, k);
    double full = c147(l, c1v, c2v, ph, p, k);
    reduceErr = std::max(reduceErr, std::fabs(red - ref) / (1e-300 + std::fabs(ref)));
    full478Diff = std::max(full478Diff, std::fabs(full - ref) / (1e-300 + std::fabs(ref)));
  }
  std::printf("ZA4 1/4/7 reduces to struct-1 (ZAqbq4=ZAqbq7=0): max rel err = %.3e\n", reduceErr);
  std::printf("ZA4 1/4/7 with real ZAqbq4/7 differs from struct-1: max rel diff = %.3e\n", full478Diff);
  // reduceErr is pure round-off: the 1/4/7 (structures off) and struct-1 kernels are independent
  // CSE/Horner straight-line programs for the SAME value, differing only by FP reassociation.
  bool ok = reduceErr < 1e-9 && full478Diff > 1e-3;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
