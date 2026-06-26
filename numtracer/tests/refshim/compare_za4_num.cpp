// Validate + TIME the NUMERIC backend for the QCD ZA4 four-gluon vertex (quark box) vs the FORM
// oracle. 3 loop angles. Physical check = the angular-integrated value; timing is ns/eval over
// random points.
#include "ZA4_kernel.hh"      // copied FormTracer oracle (refshim/)
#include "ZA4_num_kernel.hh"  // numeric backend         (gen/)
#include "shim.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

int main() {
  using Form = DiFfRG::ZA4_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZA4_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  const double Nf = 2.0;
  auto cForm = [&](double l1, double c1, double c2, double ph, double p, double k) {
    return Form::kernel(l1, c1, c2, ph, p, k, Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq);
  };
  auto cNum = [&](double l1, double c1, double c2, double ph, double p, double k) {
    return std::real(Num::kernel(l1, c1, c2, ph, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999), Uphi(0.0, 6.283185307);
  const int Np = 200000;
  std::vector<double> L(Np), C1(Np), C2(Np), PH(Np), P(Np), K(Np);
  for (int i = 0; i < Np; ++i) { L[i] = U(rng); C1[i] = Uc(rng); C2[i] = Uc(rng); PH[i] = Uphi(rng); P[i] = U(rng); K[i] = U(rng); }

  double pw = 0;
  for (int i = 0; i < Np; ++i) {
    double r = cForm(L[i], C1[i], C2[i], PH[i], P[i], K[i]);
    pw = std::max(pw, std::fabs(cNum(L[i], C1[i], C2[i], PH[i], P[i], K[i]) - r) / (1e-300 + std::fabs(r)));
  }
  auto ang3 = [&](auto fn, double l1, double p, double k) {
    const int M = 17, Mp = 24; double s = 0, h = 2.0 / (M - 1), hp = 6.283185307 / Mp;
    for (int i = 0; i < M; ++i) for (int j = 0; j < M; ++j) for (int m = 0; m < Mp; ++m) {
      double c1 = -1 + i * h, c2 = -1 + j * h, ph = m * hp;
      double w = ((i == 0 || i == M - 1) ? 0.5 : 1.0) * ((j == 0 || j == M - 1) ? 0.5 : 1.0);
      s += w * fn(l1, c1, c2, ph, p, k);
    }
    return s * h * h * hp;
  };
  double ie = 0;
  for (double l1 : {0.6, 1.3}) for (double p : {0.7, 1.8}) for (double k : {0.5, 1.4}) {
    double r = ang3(cForm, l1, p, k);
    ie = std::max(ie, std::fabs(ang3(cNum, l1, p, k) - r) / (1e-300 + std::fabs(r)));
  }
  std::printf("ZA4 numeric vs FORM:  pointwise=%.3e  3-angle-integrated=%.3e (physical check)\n", pw, ie);

  // --- timing: ns/eval over the random points ---
  auto bench = [&](auto fn, const char *nm) {
    volatile double sink = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < Np; ++i) sink = sink + fn(L[i], C1[i], C2[i], PH[i], P[i], K[i]);
    auto t1 = std::chrono::steady_clock::now();
    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / Np;
    std::printf("  %-6s %8.1f ns/eval\n", nm, ns);
    return ns;
  };
  std::printf("timing over %d random points:\n", Np);
  double nF = bench(cForm, "Form");
  double nN = bench(cNum, "num");
  std::printf("RESULT kernel=ZA4 path=num ns=%.1f relerr=%.3e  (Form=%.1f)\n", nN, ie, nF);

  bool ok = ie < 1e-8;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
