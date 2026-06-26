// Validate the NUMERIC backend for the pure-YM ZAcbc ghost-gluon vertex vs the DENSE oracle
// (entry-for-entry contraction, the trusted reference that matches FormTracer). This vertex has a
// diagram with a DISCONNECTED Lorentz network (two index-disjoint trace factors that MULTIPLY); it
// is the regression guard for the codegen bug where the numeric backend SUMMED such components
// instead of multiplying them. Numeric and Dense share the same frame, so they must agree POINTWISE.
#include "ZAcbc_num_dense_kernel.hh" // dense oracle  (gen/)
#include "ZAcbc_num_kernel.hh"       // numeric backend (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using Dense = DiFfRG::ZAcbc_num_dense_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZAcbc_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  auto cDense = [&](double l1, double c1, double c2, double p, double k) {
    return std::real(Dense::kernel(l1, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.dtZc, d.Zc, d.dtZA, d.ZA));
  };
  auto cNum = [&](double l1, double c1, double c2, double p, double k) {
    return std::real(Num::kernel(l1, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.dtZc, d.Zc, d.dtZA, d.ZA));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  const int Np = 200000;
  double pw = 0;
  for (int i = 0; i < Np; ++i) {
    double l1 = U(rng), c1 = Uc(rng), c2 = Uc(rng), p = U(rng), k = U(rng);
    double r = cDense(l1, c1, c2, p, k);
    pw = std::max(pw, std::fabs(cNum(l1, c1, c2, p, k) - r) / (1e-300 + std::fabs(r)));
  }
  auto ang2 = [&](auto kfn, double l1, double p, double k) {
    const int M = 41; double s = 0, h = 2.0 / (M - 1);
    for (int i = 0; i < M; ++i) for (int j = 0; j < M; ++j) {
      double c1 = -1.0 + i * h, c2 = -1.0 + j * h;
      double w = ((i == 0 || i == M - 1) ? 0.5 : 1.0) * ((j == 0 || j == M - 1) ? 0.5 : 1.0);
      s += w * kfn(l1, c1, c2, p, k);
    }
    return s * h * h;
  };
  double ie = 0;
  for (double l1 : {0.6, 1.3, 2.4}) for (double p : {0.7, 1.8}) for (double k : {0.5, 1.4}) {
    double r = ang2(cDense, l1, p, k), g = ang2(cNum, l1, p, k);
    ie = std::max(ie, std::fabs(g - r) / (1e-300 + std::fabs(r)));
  }
  std::printf("ZAcbc numeric vs DENSE:  pointwise=%.3e  (cos1,cos2)-integrated=%.3e\n", pw, ie);
  bool ok = (pw < 1e-9) && (ie < 1e-9);
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
