// ns/eval timing: dense 1/4/7 numeric backend kernel vs the FormTracer FORM oracle.
#include "ZA3_147_kernel.hh"     // FORM oracle (refshim/)
#include "ZA3_147_num_kernel.hh" // numeric backend (gen/)
#include "shim.hpp"
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

int main() {
  using Form = DiFfRG::ZA3_147_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZA3_147_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  const double Nf = 2.0;
  auto cForm = [&](double l, double c1, double c2, double p, double k) {
    return Form::kernel(l, c1, c2, p, k, Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                        d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq);
  };
  auto cNum = [&](double l, double c1, double c2, double p, double k) {
    return std::real(Num::kernel(l, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                                 d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  const int Np = 200000;
  std::vector<double> L(Np), C1(Np), C2(Np), P(Np), K(Np);
  for (int i = 0; i < Np; ++i) { L[i] = U(rng); C1[i] = Uc(rng); C2[i] = Uc(rng); P[i] = U(rng); K[i] = U(rng); }
  auto bench = [&](auto fn, const char *nm) {
    volatile double sink = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < Np; ++i) sink = sink + fn(L[i], C1[i], C2[i], P[i], K[i]);
    auto t1 = std::chrono::steady_clock::now();
    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / Np;
    std::printf("  %-6s %9.1f ns/eval\n", nm, ns);
    return ns;
  };
  std::printf("ZA3 1/4/7 timing over %d random points:\n", Np);
  // warm up
  bench(cForm, "Form"); bench(cNum, "num");
  double nF = bench(cForm, "Form");
  double nN = bench(cNum, "num");
  std::printf("RESULT kernel=ZA3_147  num=%.1f ns  Form=%.1f ns  ratio(num/Form)=%.2fx\n", nN, nF, nN / nF);
  return 0;
}
