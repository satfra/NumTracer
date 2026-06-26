// ZAqbq{1,4,7}_147: validate numeric == FORM and time both (ns/call), per struct.
#include "ZAqbq1_147_kernel.hh"
#include "ZAqbq4_147_kernel.hh"
#include "ZAqbq7_147_kernel.hh"
#include "ZAqbq1_147_num_kernel.hh"
#include "ZAqbq4_147_num_kernel.hh"
#include "ZAqbq7_147_num_kernel.hh"
#include "shim.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

int main() {
  DressingSet d;
  const double Nf = 2.0;
  const int Np = 200000;
  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  std::vector<std::array<double, 5>> pts(Np);
  for (auto &q : pts) q = {U(rng), Uc(rng), Uc(rng), U(rng), U(rng)};

  volatile double sink = 0;
  auto run = [&](int lab, auto formK, auto numK) {
    double pw = 0, maxabs = 0;
    for (auto &q : pts) {
      double r = formK(q[0], q[1], q[2], q[3], q[4], Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4,
                       d.ZAqbq7, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq);
      double n = std::real(numK(q[0], q[1], q[2], q[3], q[4], d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4,
                                d.ZAqbq7, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
      maxabs = std::max(maxabs, std::fabs(r));
      pw = std::max(pw, std::fabs(n - r) / (1e-300 + std::fabs(r)));
    }
    auto t0 = std::chrono::steady_clock::now();
    for (auto &q : pts)
      sink += formK(q[0], q[1], q[2], q[3], q[4], Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4,
                    d.ZAqbq7, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq);
    auto t1 = std::chrono::steady_clock::now();
    for (auto &q : pts)
      sink += std::real(numK(q[0], q[1], q[2], q[3], q[4], d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4,
                             d.ZAqbq7, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
    auto t2 = std::chrono::steady_clock::now();
    double fns = std::chrono::duration<double, std::nano>(t1 - t0).count() / Np;
    double nns = std::chrono::duration<double, std::nano>(t2 - t1).count() / Np;
    std::printf("ZAqbq%d  rel-err %.3e (|max|=%.2e)  | FORM %6.1f ns  numeric %6.1f ns  (num/FORM %.2fx)\n",
                lab, pw, maxabs, fns, nns, nns / fns);
  };

  run(1, &DiFfRG::ZAqbq1_147_kernel<DiFfRG::ShimRegulator>::kernel,
         &DiFfRG::ZAqbq1_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  run(4, &DiFfRG::ZAqbq4_147_kernel<DiFfRG::ShimRegulator>::kernel,
         &DiFfRG::ZAqbq4_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  run(7, &DiFfRG::ZAqbq7_147_kernel<DiFfRG::ShimRegulator>::kernel,
         &DiFfRG::ZAqbq7_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  return 0;
}
