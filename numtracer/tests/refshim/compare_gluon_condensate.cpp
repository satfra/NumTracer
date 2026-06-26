// DEMONSTRATION + validation — a PER-COMPONENT ADJOINT dressing (the gluon-condensation use case).
//
// The _diag kernel (gen/gen_gluon_condensate.wls) carries an SU(3)-adjoint internal line dressed
// component-by-component, diag(ZA[0],..,ZA[7]); the adjoint trace folds (via sun_value_dressed) to
// Σ_{a=0}^{7} of the per-component dressing product — a runtime adjoint-sum, with the Dirac trace
// computed ONCE (not 8 diagrams). Validation needs no dense oracle:
//   (1) COLLAPSE: feeding all 8 components the SAME function ZAflat reproduces the colour-blind
//       _plain kernel (δ^{aa}=8) exactly — Σ_a Z^2 = 8 Z^2.
//   (2) PER-COMPONENT: with distinct ZA[a] the result genuinely differs from the collapsed one,
//       i.e. the 8 colours are dressed independently (not folded to one blind constant).
#include "Gluon_condensate_diag_kernel.hh"  // per-component adjoint dressing (8 components)  (gen/)
#include "Gluon_condensate_plain_kernel.hh" // colour-blind δ^{ab} × scalar dressing          (gen/)
#include "shim.hpp"                           // DiFfRG::Fn, ShimRegulator

#include <array>
#include <cmath>
#include <cstdio>
#include <random>

namespace {
double zflat(double x) { return 1.0 / (1.0 + 0.2 * x) + 0.85; }      // the blind reference dressing
double zflatDot(double x) { return 0.50 + 0.10 * std::sin(0.3 * x); }
// 8 DISTINCT per-colour dressings (a genuinely component-dependent condensate).
double zc(double x, double s) { return 0.6 + 0.4 * std::sin(0.2 * x + s) + 0.05 * s; }
} // namespace

int main()
{
  using Diag = DiFfRG::Gluon_condensate_diag_kernel<DiFfRG::ShimRegulator>;
  using Plain = DiFfRG::Gluon_condensate_plain_kernel<DiFfRG::ShimRegulator>;

  const DiFfRG::Fn ZAflat{&zflat}, ZAflatDot{&zflatDot};
  // (1) collapsed: every adjoint component = the same blind function.
  std::array<DiFfRG::Fn, 8> ZAc, ZAcDot;
  for (auto &f : ZAc) f = ZAflat;
  for (auto &f : ZAcDot) f = ZAflatDot;
  // (2) distinct per-component dressings (a genuinely colour-dependent condensate): 8 free functions.
  std::array<DiFfRG::Fn, 8> ZAv, ZAvDot;
  ZAv[0] = {+[](double x) { return zc(x, 0.0); }};
  ZAv[1] = {+[](double x) { return zc(x, 0.5); }};
  ZAv[2] = {+[](double x) { return zc(x, 1.0); }};
  ZAv[3] = {+[](double x) { return zc(x, 1.5); }};
  ZAv[4] = {+[](double x) { return zc(x, 2.0); }};
  ZAv[5] = {+[](double x) { return zc(x, 2.5); }};
  ZAv[6] = {+[](double x) { return zc(x, 3.0); }};
  ZAv[7] = {+[](double x) { return zc(x, 3.5); }};
  for (int a = 0; a < 8; ++a) ZAvDot[a] = ZAflatDot;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);

  double collapseRel = 0.0, breakRel = 0.0, maxAbs = 0.0;
  const int N = 200000;
  for (int i = 0; i < N; ++i) {
    const double l1 = U(rng), c1 = Uc(rng), p = U(rng);
    const double diagC = static_cast<double>(Diag::kernel(l1, c1, p, ZAc, ZAcDot));   // collapsed
    const double plain = static_cast<double>(Plain::kernel(l1, c1, p, ZAflat, ZAflatDot));
    const double diagV = static_cast<double>(Diag::kernel(l1, c1, p, ZAv, ZAvDot));   // per-component
    collapseRel = std::max(collapseRel, std::fabs(diagC - plain) / (1e-300 + std::fabs(plain)));
    breakRel = std::max(breakRel, std::fabs(diagV - plain) / (1e-300 + std::fabs(plain)));
    maxAbs = std::max(maxAbs, std::fabs(diagV));
  }

  std::printf("[per-component ADJOINT dressing — gluon-condensation Z_A^a folded through SU(3)]\n");
  std::printf("  (1) all 8 components equal == colour-blind δ^{aa} kernel   max rel err = %.3e\n", collapseRel);
  std::printf("  (2) distinct per-colour dressing differs from blind        max rel = %.3e  (>0 => 8 colours independent)\n",
              breakRel);
  std::printf("  max |kernel| = %.3e\n", maxAbs);

  const bool ok = collapseRel < 1e-10 && breakRel > 1e-2 && maxAbs > 1e-6;
  std::printf("\n%s\n", ok ? "ALL TESTS PASSED" : "TESTS FAILED");
  return ok ? 0 : 1;
}
