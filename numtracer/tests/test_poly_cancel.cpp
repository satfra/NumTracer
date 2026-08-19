// Correctness gate for the inverse-atom CANCELLATION passes in numeric/mpoly.hpp:
//   divThroughMonomialAtoms — denominator is a single monomial (bare/unit-direction loop)
//   divThroughPolyAtoms     — denominator is a genuine polynomial (shifted line, k² non-monomial)
//
// These two decide which 1/k² atoms survive to become a runtime division in the emitted kernel, and
// they were the only major engine transforms with NO direct test: they appeared solely inside the
// generated gen_*.cpp, which are untracked build intermediates. A change there is therefore invisible
// to `ctest` — it shows up as a wrong or bloated kernel much later, if at all. Measured stakes:
// divThroughPolyAtoms is worth 1.42x emitted SSA on the dressed aqbq147 flows.
//
// The invariants asserted here are the ones an "improvement" to these passes actually threatens:
//   1. VALUE PRESERVATION. p and its reduction must agree numerically at random points, with each
//      surviving atom evaluated as 1/D_aid. This is the whole contract; everything else is
//      optimisation. Checked for both the fires and the does-not-fire cases.
//   2. THE CANCELLATION ACTUALLY FIRES on an exactly-divisible numerator (atom gone, value = Q).
//      Without this a pass that returns its input unchanged would satisfy invariant 1 forever.
//   3. NO TERM-COUNT GROWTH. These are simplification passes; a reduction that returns more
//      monomials than it consumed is a pessimisation. This is not hypothetical — a partial-fraction
//      relaxation of divThroughPolyAtoms (N/D -> Q + R/D) was implemented and measured 2026-08-08 at
//      1.5-2.5x MORE emitted SSA on 8/8 flows, precisely because it splits one group into two. It
//      passed value preservation the whole time. Term count is the invariant that catches it.
//   4. IDEMPOTENCE. Both passes are documented as running to a fixed point.
#include "numtracer/numeric/env.hpp"
#include "numtracer/numeric/mpoly.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace nm = numtracer::numeric;
using numtracer::Cx;

static int g_fails = 0;

static void check(bool ok, const char *what)
{
  if (!ok) {
    std::printf("  FAIL  %s\n", what);
    ++g_fails;
  }
}

/// Count monomials still carrying at least one inverse atom.
static std::size_t termsWithAtoms(const nm::MPoly &p)
{
  std::size_t n = 0;
  for (const auto &[m, c] : p.terms) {
    (void)c;
    if (!m.atoms.empty()) ++n;
  }
  return n;
}

/// Evaluate `p` at `x`, with atom `aid` standing for 1/atomDen[aid](x) — the semantics the passes
/// must preserve. Returns false if any denominator is (near) zero at this point.
static bool evalAt(const nm::MPoly &p, const std::vector<double> &x, const std::vector<nm::MPoly> &atomDen, Cx &out)
{
  std::vector<double> inv(atomDen.size(), 0.0);
  for (std::size_t a = 0; a < atomDen.size(); ++a) {
    const Cx d = nm::eval(atomDen[a], x, {});
    if (std::fabs(d.re) < 1e-3 || std::fabs(d.im) > 1e-12) return false;
    inv[a] = 1.0 / d.re;
  }
  out = nm::eval(p, x, inv);
  return true;
}

static double relerr(Cx a, Cx b)
{
  const double s = std::max(1.0, std::max(std::fabs(a.re), std::fabs(a.im)));
  return std::max(std::fabs(a.re - b.re), std::fabs(a.im - b.im)) / s;
}

/// Assert value preservation of `red` vs `orig` over many random points.
static void checkSameValue(const char *what, const nm::MPoly &orig, const nm::MPoly &red,
                           const std::vector<nm::MPoly> &atomDen, int nsym, std::mt19937 &rng)
{
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  double worst = 0.0;
  int pts = 0;
  for (int trial = 0; trial < 400 && pts < 60; ++trial) {
    std::vector<double> x(static_cast<std::size_t>(nsym));
    for (double &v : x)
      v = U(rng);
    Cx a{}, b{};
    if (!evalAt(orig, x, atomDen, a)) continue;
    if (!evalAt(red, x, atomDen, b)) continue;
    ++pts;
    worst = std::max(worst, relerr(a, b));
  }
  check(pts >= 20, "enough well-conditioned sample points");
  check(worst < 1e-9, what);
  if (worst >= 1e-9) std::printf("        worst relative deviation %.3e over %d points\n", worst, pts);
}

int main()
{
  std::mt19937 rng(20260808);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  const int nsym = 4;
  nm::LorentzEnv env(nsym);

  auto x0 = env.var(0), x1 = env.var(1), x2 = env.var(2), x3 = env.var(3);

  // A multi-term, atom-free denominator: the shifted-line k² case divThroughPolyAtoms exists for.
  const nm::MPoly D = x0 * x0 + x1 * x1 + env.constant(Cx{2.0, 0.0}) * x0 * x1 + env.constant(Cx{3.0, 0.0});
  // A single-monomial denominator: the bare-loop case divThroughMonomialAtoms exists for.
  const nm::MPoly Dm = x2 * x2;

  // ---- 1) divThroughPolyAtoms FIRES on an exactly divisible numerator -------------------------
  {
    std::printf("== divThroughPolyAtoms: exact division cancels the atom ==\n");
    const std::vector<nm::MPoly> aden{D};
    const nm::MPoly Q = x0 * x1 + x3 * x3 + env.constant(Cx{-1.5, 0.0});
    const nm::MPoly p = (Q * D) * env.atom(0); // (Q·D)·(1/D)
    check(termsWithAtoms(p) == p.terms.size() && !p.terms.empty(), "input carries the atom on every term");

    const nm::MPoly r = nm::divThroughPolyAtoms(p, aden);
    check(termsWithAtoms(r) == 0, "atom fully cancelled on exact division");
    check(r.terms.size() == Q.terms.size(), "reduced to the quotient's term count");
    checkSameValue("exact-division value preserved", p, r, aden, nsym, rng);
    check(r.terms.size() <= p.terms.size(), "no term-count growth");

    const nm::MPoly r2 = nm::divThroughPolyAtoms(r, aden);
    check(r2.terms.size() == r.terms.size(), "idempotent");
  }

  // ---- 2) divThroughPolyAtoms on a NON-divisible numerator ------------------------------------
  // The atom must survive and the value must be untouched. This is where a partial-fraction
  // relaxation would silently change the shape; invariant 3 is what rejects it.
  {
    std::printf("== divThroughPolyAtoms: non-divisible numerator is value-preserved ==\n");
    const std::vector<nm::MPoly> aden{D};
    const nm::MPoly N = x0 * x2 + x3 + env.constant(Cx{0.7, 0.0}); // deliberately not a multiple of D
    const nm::MPoly p = N * env.atom(0);

    const nm::MPoly r = nm::divThroughPolyAtoms(p, aden);
    checkSameValue("non-divisible value preserved", p, r, aden, nsym, rng);
    check(r.terms.size() <= p.terms.size(), "no term-count growth on a failed division");
    check(termsWithAtoms(r) > 0, "the atom survives (it does not divide out)");
  }

  // ---- 3) divThroughMonomialAtoms: term-by-term cancellation ----------------------------------
  {
    std::printf("== divThroughMonomialAtoms: monomial denominator ==\n");
    const std::vector<nm::MPoly> aden{Dm};
    // x2^2 · (1/x2^2) must cancel exactly; the x3 term cannot and must keep its atom.
    const nm::MPoly p = (x2 * x2 * x0) * env.atom(0) + (x3 * x1) * env.atom(0);
    const nm::MPoly r = nm::divThroughMonomialAtoms(p, aden);
    checkSameValue("monomial-cancellation value preserved", p, r, aden, nsym, rng);
    check(r.terms.size() <= p.terms.size(), "no term-count growth");
    check(termsWithAtoms(r) < termsWithAtoms(p), "at least one atom cancelled");

    const nm::MPoly r2 = nm::divThroughMonomialAtoms(r, aden);
    check(r2.terms.size() == r.terms.size(), "idempotent");
  }

  // ---- 4) randomised value preservation over mixed inputs -------------------------------------
  // Both passes composed, on numerators that are sometimes divisible and usually not.
  {
    std::printf("== randomised: composed passes preserve value ==\n");
    const std::vector<nm::MPoly> aden{D, Dm};
    int cases = 0;
    for (int it = 0; it < 40; ++it) {
      nm::MPoly N = env.constant(Cx{U(rng), 0.0});
      for (int k = 0; k < 3; ++k) {
        const int i = static_cast<int>((rng() % 4));
        const int j = static_cast<int>((rng() % 4));
        N = N + env.constant(Cx{U(rng), 0.0}) * env.var(i) * env.var(j);
      }
      if (it % 3 == 0) N = N * D; // make a third of them exactly divisible
      const int aid = static_cast<int>(it % 2);
      const nm::MPoly p = N * env.atom(aid);

      nm::MPoly r = nm::divThroughMonomialAtoms(p, aden);
      r = nm::divThroughPolyAtoms(r, aden);
      checkSameValue("composed value preserved", p, r, aden, nsym, rng);
      check(r.terms.size() <= p.terms.size(), "no term-count growth (composed)");
      ++cases;
    }
    check(cases == 40, "all randomised cases ran");
  }

  if (g_fails) {
    std::printf("\ntest_poly_cancel: %d FAILURE(S)\n", g_fails);
    return 1;
  }
  std::printf("\ntest_poly_cancel: all checks passed\n");
  return 0;
}
