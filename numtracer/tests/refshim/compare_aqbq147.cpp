// Validate the dense 1/4/7 NUMERIC backend for the QCD quark-gluon vertex flow (A qb q 3-point),
// projected onto EACH external tensor structure (1, 4, 7) of the AqbqDirect147 basis, against the
// copied FormTracer oracles. Like the ZAqbq1 vertex flow these match the oracle POINTWISE (no odd
// loop-routing term that integrates to zero). The numeric kernels carry the fixed-frame round-off
// floor, so the pointwise tolerance is 1e-8 (as for the other numeric kernels).
#include "ZAqbq1_147_kernel.hh"     // FormTracer oracle, proj struct 1 (refshim/)
#include "ZAqbq4_147_kernel.hh"     // FormTracer oracle, proj struct 4 (refshim/)
#include "ZAqbq7_147_kernel.hh"     // FormTracer oracle, proj struct 7 (refshim/)
#include "ZAqbq1_147_num_kernel.hh" // numeric backend, proj struct 1 (gen/)
#include "ZAqbq4_147_num_kernel.hh" // numeric backend, proj struct 4 (gen/)
#include "ZAqbq7_147_num_kernel.hh" // numeric backend, proj struct 7 (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main() {
  DressingSet d;
  const double Nf = 2.0;

  // each structure: {label, FORM kernel, numeric kernel}
  auto run = [&](int lab, auto formK, auto numK) {
    auto cForm = [&](double l1, double c1, double c2, double p, double k) {
      return formK(l1, c1, c2, p, k, Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                   d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq);
    };
    auto cNum = [&](double l1, double c1, double c2, double p, double k) {
      return std::real(numK(l1, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7,
                            d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
    };
    std::mt19937_64 rng(31337);
    std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
    const int Np = 200000;
    double pw = 0, maxabs = 0;
    for (int i = 0; i < Np; ++i) {
      double l1 = U(rng), c1 = Uc(rng), c2 = Uc(rng), p = U(rng), k = U(rng);
      double r = cForm(l1, c1, c2, p, k);
      maxabs = std::max(maxabs, std::fabs(r));
      pw = std::max(pw, std::fabs(cNum(l1, c1, c2, p, k) - r) / (1e-300 + std::fabs(r)));
    }
    bool ok = pw < 1e-8;
    std::printf("  ZAqbq%d (147) numeric vs FORM:  pointwise rel err = %.3e  (|max|=%.2e)  %s\n",
                lab, pw, maxabs, ok ? "ok" : "FAIL");
    return ok;
  };

  std::printf("QCD quark-gluon vertex flow, full {1,4,7} basis, projected per structure:\n");
  bool ok = true;
  ok &= run(1, &DiFfRG::ZAqbq1_147_kernel<DiFfRG::ShimRegulator>::kernel,
               &DiFfRG::ZAqbq1_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  ok &= run(4, &DiFfRG::ZAqbq4_147_kernel<DiFfRG::ShimRegulator>::kernel,
               &DiFfRG::ZAqbq4_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  ok &= run(7, &DiFfRG::ZAqbq7_147_kernel<DiFfRG::ShimRegulator>::kernel,
               &DiFfRG::ZAqbq7_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  std::printf(ok ? "\nALL TESTS PASSED\n" : "\nTESTS FAILED\n");
  return ok ? 0 : 1;
}
