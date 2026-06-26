// ns/eval timing of the DENSE oracle (entry-for-entry brute force) on the QCD ZA gluon
// self-energy — fast enough (~ms/eval) to measure the effect of changes to the dense fold
// (dense/dtensor.hpp). Also cross-checks the dense oracle against the numeric backend.
//
//   cmake -S . -B build -DNUMTRACER_BUILD_ZA3_147_DENSE=ON
//   cmake --build build --target bench_za_dense && NT_DENSE_NP=20000 ./build/bench_za_dense
#include "ZA_num_dense_kernel.hh" // dense oracle    (gen/)
#include "ZA_num_kernel.hh"       // numeric backend (gen/)
#include "shim.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

int main() {
  using Dense = DiFfRG::ZA_num_dense_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZA_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  auto cDense = [&](double l1, double c1, double p, double k) {
    return std::real(Dense::kernel(l1, c1, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  auto cNum = [&](double l1, double c1, double p, double k) {
    return std::real(Num::kernel(l1, c1, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  int Np = 5000;
  if (const char *e = std::getenv("NT_DENSE_NP")) Np = std::max(1, std::atoi(e));
  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  std::vector<double> L(Np), C1(Np), P(Np), K(Np);
  for (int i = 0; i < Np; ++i) { L[i] = U(rng); C1[i] = Uc(rng); P[i] = U(rng); K[i] = U(rng); }

  double maxrel = 0.0;
  for (int i = 0; i < std::min(Np, 2000); ++i) {
    double rel = std::fabs(cDense(L[i], C1[i], P[i], K[i]) - cNum(L[i], C1[i], P[i], K[i])) /
                 (1e-30 + std::fabs(cNum(L[i], C1[i], P[i], K[i])));
    maxrel = std::max(maxrel, rel);
  }
  std::printf("ZA dense-vs-numeric: max rel err = %.3e\n", maxrel);

  auto bench = [&](auto fn, const char *nm) {
    volatile double sink = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < Np; ++i) sink = sink + fn(L[i], C1[i], P[i], K[i]);
    auto t1 = std::chrono::steady_clock::now();
    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / Np;
    std::printf("  %-6s %12.1f ns/eval\n", nm, ns);
    return ns;
  };
  std::printf("ZA dense timing over %d random points:\n", Np);
  bench(cDense, "dense"); // warm up
  double nD = bench(cDense, "dense");
  double nN = bench(cNum, "num");
  std::printf("RESULT kernel=ZA_dense  dense=%.1f ns  num=%.1f ns  ratio(dense/num)=%.2fx\n", nD, nN, nD / nN);
  return maxrel < 1e-9 ? 0 : 1;
}
