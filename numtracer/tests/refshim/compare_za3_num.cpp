// Validate the NUMERIC backend for the QCD ZA3 three-gluon vertex (quark triangle) vs the copied
// FormTracer oracle. 2 loop angles. As for inv, pointwise can differ by an odd loop-routing term
// that integrates to zero — the physical check is the (cos1,cos2)-integrated value.
#include "ZA3_kernel.hh"     // copied FormTracer oracle (refshim/)
#include "ZA3_num_kernel.hh" // numeric backend         (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using Form = DiFfRG::ZA3_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZA3_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  const double Nf = 2.0;
  auto cForm = [&](double l1, double c1, double c2, double p, double k) {
    return Form::kernel(l1, c1, c2, p, k, Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq);
  };
  auto cNum = [&](double l1, double c1, double c2, double p, double k) {
    return std::real(Num::kernel(l1, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  const int Np = 200000;
  double pw = 0;
  for (int i = 0; i < Np; ++i) {
    double l1 = U(rng), c1 = Uc(rng), c2 = Uc(rng), p = U(rng), k = U(rng);
    double r = cForm(l1, c1, c2, p, k);
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
    double r = ang2(cForm, l1, p, k), g = ang2(cNum, l1, p, k);
    ie = std::max(ie, std::fabs(g - r) / (1e-300 + std::fabs(r)));
  }
  std::printf("ZA3 numeric vs FORM:  pointwise=%.3e  (cos1,cos2)-integrated=%.3e (physical check)\n", pw, ie);
  bool ok = ie < 1e-9;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
