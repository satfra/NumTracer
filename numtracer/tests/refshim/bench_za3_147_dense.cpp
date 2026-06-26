// Stack-overflow reproducer + ns/eval timing harness for the DENSE backend, on the largest
// committed dense kernel (ZA3 1/4/7). On the unfixed backend this SIGSEGVs (a single huge
// contract_all intermediate is a stack std::array); after the heap-backed fold it runs to
// completion. Cross-checks the dense oracle against the numeric backend over many random points.
//
//   cmake -S . -B build -DNUMTRACER_BUILD_ZA3_147_DENSE=ON
//   cmake --build build --target bench_za3_147_dense && ./build/bench_za3_147_dense
#include "ZA3_147_num_dense_kernel.hh" // dense oracle (gen/)
#include "ZA3_147_num_kernel.hh"       // numeric backend (gen/)
#include "shim.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

int main() {
  using Dense = DiFfRG::ZA3_147_num_dense_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZA3_147_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  auto cDense = [&](double l, double c1, double c2, double p, double k) {
    return std::real(Dense::kernel(l, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                                   d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  auto cNum = [&](double l, double c1, double c2, double p, double k) {
    return std::real(Num::kernel(l, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                                 d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  // The dense ZA3 1/4/7 oracle is ~30-60s PER eval (it sums every index combination of every
  // diagram of a huge loop), so this driver exists to prove the heap-backed fold no longer
  // overflows the stack and still matches the numeric backend — not as a timing vehicle (use
  // bench_za_dense for that). A handful of points suffices; NT_DENSE_NP overrides.
  int Np = 3;
  if (const char *e = std::getenv("NT_DENSE_NP")) Np = std::max(1, std::atoi(e));
  std::vector<double> L(Np), C1(Np), C2(Np), P(Np), K(Np);
  for (int i = 0; i < Np; ++i) { L[i] = U(rng); C1[i] = Uc(rng); C2[i] = Uc(rng); P[i] = U(rng); K[i] = U(rng); }

  // Correctness: the dense oracle must match the numeric backend pointwise.
  double maxrel = 0.0;
  int checks = std::min(Np, 20);
  for (int i = 0; i < checks; ++i) {
    double dv = cDense(L[i], C1[i], C2[i], P[i], K[i]);
    double nv = cNum(L[i], C1[i], C2[i], P[i], K[i]);
    double rel = std::abs(dv - nv) / (std::abs(nv) + 1e-30);
    if (rel > maxrel) maxrel = rel;
  }
  std::printf("ZA3 1/4/7 dense-vs-numeric over %d points: max rel err = %.3e\n", checks, maxrel);

  auto bench = [&](auto fn, const char *nm) {
    volatile double sink = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < Np; ++i) sink = sink + fn(L[i], C1[i], C2[i], P[i], K[i]);
    auto t1 = std::chrono::steady_clock::now();
    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / Np;
    std::printf("  %-6s %12.1f ns/eval\n", nm, ns);
    return ns;
  };
  std::printf("ZA3 1/4/7 dense timing over %d random points:\n", Np);
  double nD = bench(cDense, "dense");
  double nN = bench(cNum, "num");
  std::printf("RESULT kernel=ZA3_147_dense  dense=%.1f ns  num=%.1f ns  ratio(dense/num)=%.2fx\n", nD, nN, nD / nN);
  return maxrel < 1e-9 ? 0 : 1;
}
