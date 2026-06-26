// Validate the NUMERIC (matrix-product) backend for the QCD ZA gluon self-energy against the Dense
// oracle (the entry-for-entry baseline, itself validated vs FORM). Pointwise + cos1-integrated.
#include "ZA_num_dense_kernel.hh" // dense oracle   (gen/)
#include "ZA_num_kernel.hh"       // numeric backend (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using Dense = DiFfRG::ZA_num_dense_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZA_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  auto callDense = [&](double l1, double c1, double p, double k) {
    return std::real(Dense::kernel(l1, c1, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  auto callNum = [&](double l1, double c1, double p, double k) {
    return std::real(Num::kernel(l1, c1, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  const int Np = 200000;
  double pw = 0;
  for (int i = 0; i < Np; ++i) {
    double l1 = U(rng), c1 = Uc(rng), p = U(rng), k = U(rng);
    double r = callDense(l1, c1, p, k), g = callNum(l1, c1, p, k);
    pw = std::max(pw, std::fabs(g - r) / (1e-300 + std::fabs(r)));
  }
  auto ang1 = [&](auto kfn, double l1, double p, double k) {
    const int M = 321; double s = 0, h = 2.0 / (M - 1);
    for (int i = 0; i < M; ++i) { double c1 = -1.0 + i * h, w = (i == 0 || i == M - 1) ? 0.5 : 1.0; s += w * kfn(l1, c1, p, k); }
    return s * h;
  };
  double ie = 0;
  for (double l1 : {0.6, 1.3, 2.4}) for (double p : {0.7, 1.8}) for (double k : {0.5, 1.4}) {
    double r = ang1(callDense, l1, p, k), g = ang1(callNum, l1, p, k);
    ie = std::max(ie, std::fabs(g - r) / (1e-300 + std::fabs(r)));
  }
  std::printf("ZA numeric vs dense:  pointwise=%.3e  cos1-integrated=%.3e\n", pw, ie);
  bool ok = pw < 1e-10 && ie < 1e-10;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
