// Symbolic dressing collection gate: the Zq quark self-energy generated with DressingCollection ->
// True (the quark propagator numerator Mq·δ + Zq·γ·l kept EAGER as one ntDressedNum slot, folded to a
// single DPoly trace) must reproduce the committed DISTRIBUTED Zq kernel (which is itself validated vs
// the FORM reference by flow_zq_num) to round-off. This exercises the full dressed path — front-end
// ntDressedNum rewrite (with common colour/flavour/denominator factored out), the DPoly generator
// branch, the kind-2 `dress` env leaves, and the dressed contraction — on a real flow.
#include "Zq_collect_kernel.hh"
#include "Zq_num_kernel.hh"
#include "shim.hpp"
#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using Dist = DiFfRG::Zq_num_kernel;      // distributed (FORM-validated)
  using Coll = DiFfRG::Zq_collect_kernel;  // symbolic dressing collection
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
  // PHYSICAL CHECK = the angular-integrated value. The pointwise integrand is not a physical
  // invariant: the loop momentum is a free reparametrisation, and FunKit's FRoute picks it from the
  // leg order it is handed — so a change of backend or derivation can hand us an integrand that
  // differs pointwise by O(1) while integrating to exactly the same thing. Integrate the one loop
  // angle with its sqrt(1-cos1^2) polar Jacobian, at fixed |l1|.
  auto ang1 = [&](auto kfn, double l1, double p, double k) {
    const int M = 65; double s = 0, h = 2.0 / (M - 1);
    for (int i = 0; i < M; ++i) {
      double c1 = -1.0 + i * h;
      double w = ((i == 0 || i == M - 1) ? 0.5 : 1.0);
      s += w * std::sqrt(std::max(0.0, 1.0 - c1 * c1)) * kfn(l1, c1, p, k);
    }
    return s * h;
  };
  double ie = 0;
  for (double l1 : {0.4, 0.9, 1.7}) for (double p : {0.7, 1.8}) for (double k : {0.5, 1.4}) {
    double r = ang1(cDist, l1, p, k);
    ie = std::max(ie, std::fabs(ang1(cColl, l1, p, k) - r) / (1e-300 + std::fabs(r)));
  }
  std::printf("Zq collected vs distributed:  pointwise=%.3e  cos1-integrated=%.3e (physical check)\n", pw, ie);
  bool ok = ie < 1e-8; // numeric (evaluation-order round-off) tolerance, as in compare_zq_num
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
