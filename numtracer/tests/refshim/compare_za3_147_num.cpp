// Cross-check the dense 1/4/7 numeric kernel against the FORM-validated struct-1 numeric kernel.
// With ZAqbq4 = ZAqbq7 = 0 the 1/4/7 vertex collapses to structure 1, so the 1/4/7 kernel must
// reduce EXACTLY (round-off) to ZA3_num_kernel (which compare_za3_num validates vs FORM). With the
// real ZAqbq4/ZAqbq7 it must DIFFER (structures 4/7 genuinely contribute). Both kernels are the
// compact numeric backend (fast), so no Dense oracle / huge stack is needed.
#include "ZA3_147_num_kernel.hh" // dense 1/4/7 numeric kernel (gen/)
#include "ZA3_num_kernel.hh"     // struct-1 numeric kernel, validated vs FORM (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using K147 = DiFfRG::ZA3_147_num_kernel<DiFfRG::ShimRegulator>;
  using K1 = DiFfRG::ZA3_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  DiFfRG::Fn zero{+[](double) { return 0.0; }};

  // 1/4/7 kernel with structures 4,7 turned off:
  auto c147_0 = [&](double l, double c1, double c2, double p, double k) {
    return std::real(K147::kernel(l, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, zero, zero,
                                  d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  // 1/4/7 kernel with the real structure-4,7 dressings:
  auto c147 = [&](double l, double c1, double c2, double p, double k) {
    return std::real(K147::kernel(l, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                                  d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  // struct-1 kernel (the FORM-validated reference):
  auto c1 = [&](double l, double cc1, double cc2, double p, double k) {
    return std::real(K1::kernel(l, cc1, cc2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1,
                                d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(20240618);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  const int Np = 200000;
  double reduceErr = 0, full478Diff = 0;
  for (int i = 0; i < Np; ++i) {
    double l = U(rng), c1v = Uc(rng), c2v = Uc(rng), p = U(rng), k = U(rng);
    double ref = c1(l, c1v, c2v, p, k);
    double red = c147_0(l, c1v, c2v, p, k);
    double full = c147(l, c1v, c2v, p, k);
    reduceErr = std::max(reduceErr, std::fabs(red - ref) / (1e-300 + std::fabs(ref)));
    full478Diff = std::max(full478Diff, std::fabs(full - ref) / (1e-300 + std::fabs(ref)));
  }
  std::printf("1/4/7 reduces to struct-1 (ZAqbq4=ZAqbq7=0): max rel err = %.3e\n", reduceErr);
  std::printf("1/4/7 with real ZAqbq4/7 differs from struct-1: max rel diff = %.3e\n", full478Diff);
  // reduceErr is pure round-off: the 1/4/7 (structures off) and struct-1 kernels are independent
  // CSE/Horner straight-line programs for the SAME value, so they differ only by FP reassociation
  // (~1e-11). A genuine struct-1 discrepancy would be O(1); full478Diff confirms 4/7 are present.
  bool ok = reduceErr < 1e-9 && full478Diff > 1e-3;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
