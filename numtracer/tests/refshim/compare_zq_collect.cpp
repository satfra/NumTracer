// Symbolic dressing collection gate: the Zq quark self-energy generated with DressingCollection ->
// True (the quark propagator numerator Mq·δ + Zq·γ·l kept EAGER as one ntDressedNum slot, folded to a
// single DPoly trace) must reproduce the committed DISTRIBUTED Zq kernel (which is itself validated vs
// the FORM reference by flow_zq_num) to round-off. This exercises the full dressed path — front-end
// ntDressedNum rewrite (with common colour/flavour/denominator factored out), the DPoly generator
// branch, the kind-2 `dress` env leaves, and numeric_value_dressed_netval — on a real flow.
#include "Zq_collect_kernel.hh"
#include "Zq_num_kernel.hh"
#include "shim.hpp"
#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using Dist = DiFfRG::Zq_num_kernel<DiFfRG::ShimRegulator>;      // distributed (FORM-validated)
  using Coll = DiFfRG::Zq_collect_kernel<DiFfRG::ShimRegulator>;  // symbolic dressing collection
  DressingSet d;
  auto cDist = [&](double l1, double c1, double p, double k) {
    return std::real(Dist::kernel(l1, c1, p, k, d.ZA, d.Zq, d.Mq, d.ZAqbq1, d.dtZA, d.dtZq));
  };
  auto cColl = [&](double l1, double c1, double p, double k) {
    return std::real(Coll::kernel(l1, c1, p, k, d.ZA, d.Zq, d.Mq, d.ZAqbq1, d.dtZA, d.dtZq));
  };
  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  double pw = 0;
  for (int i = 0; i < 200000; ++i) {
    double l1 = U(rng), c1 = Uc(rng), p = U(rng), k = U(rng);
    double r = cDist(l1, c1, p, k);
    pw = std::max(pw, std::fabs(cColl(l1, c1, p, k) - r) / (1e-300 + std::fabs(r)));
  }
  std::printf("Zq collected vs distributed:  pointwise=%.3e\n", pw);
  bool ok = pw < 1e-8; // numeric (evaluation-order round-off) tolerance, as in compare_zq_num
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
