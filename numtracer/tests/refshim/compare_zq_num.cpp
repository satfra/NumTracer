// Validate the NUMERIC backend for the Zq quark wave-function flow vs the Dense oracle (1 angle).
#include "Zq_num_dense_kernel.hh"
#include "Zq_num_kernel.hh"
#include "shim.hpp"
#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using Dense = DiFfRG::Zq_num_dense_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::Zq_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  auto cD = [&](double l1, double c1, double p, double k) {
    return std::real(Dense::kernel(l1, c1, p, k, d.ZA, d.Zq, d.Mq, d.ZAqbq1, d.dtZA, d.dtZq));
  };
  auto cN = [&](double l1, double c1, double p, double k) {
    return std::real(Num::kernel(l1, c1, p, k, d.ZA, d.Zq, d.Mq, d.ZAqbq1, d.dtZA, d.dtZq));
  };
  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  double pw = 0;
  for (int i = 0; i < 200000; ++i) {
    double l1 = U(rng), c1 = Uc(rng), p = U(rng), k = U(rng);
    double r = cD(l1, c1, p, k);
    pw = std::max(pw, std::fabs(cN(l1, c1, p, k) - r) / (1e-300 + std::fabs(r)));
  }
  auto ang = [&](auto fn, double l1, double p, double k) {
    const int M = 321; double s = 0, h = 2.0 / (M - 1);
    for (int i = 0; i < M; ++i) { double c1 = -1 + i * h, w = (i == 0 || i == M - 1) ? 0.5 : 1.0; s += w * fn(l1, c1, p, k); }
    return s * h;
  };
  double ie = 0;
  for (double l1 : {0.6, 1.3, 2.4}) for (double p : {0.7, 1.8}) for (double k : {0.5, 1.4}) {
    double r = ang(cD, l1, p, k);
    ie = std::max(ie, std::fabs(ang(cN, l1, p, k) - r) / (1e-300 + std::fabs(r)));
  }
  std::printf("Zq numeric vs dense:  pointwise=%.3e  cos1-integrated=%.3e\n", pw, ie);
  // Cross-compiler round-off bound: this compares the numeric kernel against the dense oracle in
  // one binary, so the gate measures FP evaluation-order round-off, not physics. Legal FMA/
  // reassociation differs by compiler (GCC ~6e-11, clang ~2e-10), so use the documented numeric
  // tolerance 1e-8 rather than a gratuitously tight bound that one compiler happens to clear.
  bool ok = pw < 1e-8 && ie < 1e-8;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
