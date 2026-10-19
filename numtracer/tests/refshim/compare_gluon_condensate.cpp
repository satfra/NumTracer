// DEMONSTRATION + validation — a per-component ADJOINT dressing (the gluon-condensation use case),
// dressing only the Cartan directions (λ3, λ8) of SU(3) and DROPPING the other 6 colours.
//
// Three kernels (gen/gen_gluon_condensate.wls) share the same physics (a scalar sigma two-point from
// a loop carrying a Dirac index and an SU(3)-adjoint index):
//   * _plain : colour-blind δ^{ab} × a single scalar dressing ZAflat (folds δ^{aa}=8).
//   * _diag  : ntSUNDiagAdj[3,a,b,{7->ZA3, 8->ZA8, Default->ZAoff}, s] — λ3 & λ8 named, the other 6
//              collapsed into a shared ZAoff (all 8 components covered).
//   * _cartan: ntSUNDiagAdj[3,a,b,{7->ZA3, 8->ZA8}, s] — only λ3 & λ8 (components 7 & 8) survive, the other 6 DROP.
// The adjoint sum folds (via sun_value_dressed) to a runtime token, with the Dirac trace computed
// ONCE (not one diagram per component). Validation needs no dense oracle:
//   (1) COLLAPSE : feeding _diag one function for ZA3=ZA8=ZAoff reproduces the blind _plain kernel
//                  exactly (Σ over 8 = 8 Z^2) — the named+Default sum covers the full group.
//   (2) DROP     : _cartan (only λ3 & λ8, engine drops the rest) equals _diag with ZAoff≡0 — the two
//                  distinct engine paths (componentDressing = -1 drop vs a zeroed Default dressing) agree.
//   (3) BREAK    : with genuinely distinct ZA3/ZA8, _cartan differs from the blind kernel, i.e. the
//                  Cartan colours are dressed independently (not folded to one blind constant).
#include "Gluon_condensate_cartan_kernel.hh" // Cartan-only (λ3 & λ8), rest dropped            (gen/)
#include "Gluon_condensate_diag_kernel.hh"   // λ3 & λ8 named + Default ZAoff for the other 6 (gen/)
#include "Gluon_condensate_plain_kernel.hh"  // colour-blind δ^{ab} × scalar dressing          (gen/)
#include "shim.hpp"                           // DiFfRG::Fn, ShimRegulator

#include <cmath>
#include <cstdio>
#include <random>

namespace {
double zflat(double x) { return 1.0 / (1.0 + 0.2 * x) + 0.85; }      // the blind reference dressing
double zflatDot(double x) { return 0.50 + 0.10 * std::sin(0.3 * x); }
// distinct Cartan dressings (a genuinely component-dependent condensate along λ3, λ8).
double z3(double x) { return 0.6 + 0.4 * std::sin(0.2 * x); }
double z8(double x) { return 0.9 - 0.3 * std::cos(0.15 * x); }
double z3dot(double x) { return 0.45 + 0.12 * std::sin(0.35 * x); }
double z8dot(double x) { return 0.55 - 0.08 * std::cos(0.28 * x); }
double zero(double) { return 0.0; }
} // namespace

int main()
{
  using Cartan = DiFfRG::Gluon_condensate_cartan_kernel;
  using Diag = DiFfRG::Gluon_condensate_diag_kernel;
  using Plain = DiFfRG::Gluon_condensate_plain_kernel;

  const DiFfRG::Fn ZAflat{&zflat}, ZAflatDot{&zflatDot};
  const DiFfRG::Fn ZA3{&z3}, ZA8{&z8}, ZA3dot{&z3dot}, ZA8dot{&z8dot};
  const DiFfRG::Fn ZERO{&zero};

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);

  double collapseRel = 0.0, dropRel = 0.0, breakRel = 0.0, maxAbs = 0.0;
  const int N = 200000;
  for (int i = 0; i < N; ++i) {
    const double l1 = U(rng), c1 = Uc(rng), p = U(rng);
    // (1) collapse: ZA3 = ZA8 = ZAoff = the blind function ⇒ _diag == _plain (Σ over 8 = 8 Z^2).
    const double diagFlat = static_cast<double>(
        Diag::kernel(l1, c1, p, ZAflat, ZAflat, ZAflat, ZAflatDot, ZAflatDot, ZAflatDot));
    const double plain = static_cast<double>(Plain::kernel(l1, c1, p, ZAflat, ZAflatDot));
    // (2) drop: _cartan (only λ3 & λ8) == _diag with the Default ZAoff ≡ 0.
    const double cartan = static_cast<double>(Cartan::kernel(l1, c1, p, ZA3, ZA8, ZA3dot, ZA8dot));
    const double diagOffZero = static_cast<double>(
        Diag::kernel(l1, c1, p, ZA3, ZA8, ZERO, ZA3dot, ZA8dot, ZERO));
    collapseRel = std::max(collapseRel, std::fabs(diagFlat - plain) / (1e-300 + std::fabs(plain)));
    dropRel = std::max(dropRel, std::fabs(cartan - diagOffZero) / (1e-300 + std::fabs(diagOffZero)));
    breakRel = std::max(breakRel, std::fabs(cartan - plain) / (1e-300 + std::fabs(plain)));
    maxAbs = std::max(maxAbs, std::fabs(cartan));
  }

  std::printf("[Cartan-dressed ADJOINT — gluon condensation along λ3, λ8 folded through SU(3)]\n");
  std::printf("  (1) ZA3=ZA8=ZAoff == colour-blind δ^{aa} kernel        max rel err = %.3e\n", collapseRel);
  std::printf("  (2) drop (cartan) == diag with ZAoff≡0                  max rel err = %.3e\n", dropRel);
  std::printf("  (3) distinct Cartan dressing differs from blind         max rel = %.3e  (>0 => independent)\n",
              breakRel);
  std::printf("  max |kernel| = %.3e\n", maxAbs);

  const bool ok = collapseRel < 1e-10 && dropRel < 1e-10 && breakRel > 1e-2 && maxAbs > 1e-6;
  std::printf("\n%s\n", ok ? "ALL TESTS PASSED" : "TESTS FAILED");
  return ok ? 0 : 1;
}
