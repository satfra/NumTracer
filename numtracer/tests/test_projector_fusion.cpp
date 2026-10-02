/// @file test_projector_fusion.cpp
/// @brief Unit tests for the projector-algebra fusion (`ndetail::fuse_projectors`): two projectors on
///        the SAME momentum contracting over a shared dummy index collapse by the algebra —
///        idempotency `P·P→P` (one shared index) / `tr P` (both shared: T→3,L→1,E→1,M→2), and
///        orthogonality `L⟂{T,E,M}`, `E⟂M` → the term vanishes. `T·E`/`T·M` are left untouched.
///
/// Two layers:
///   (A) STRUCTURAL — call fuse_projectors on a hand-built element list and assert the folded list +
///       coeff (deterministic; pins the fold logic, incl. all symmetric index-slot variants, chains,
///       traces, orthogonality, and the no-op negatives).
///   (B) VALUE — contract a redundant net vs its hand-simplified form through the contraction and check
///       they agree at random points (fusion is value-preserving).
#include <cstdio>
#include <random>
#include <vector>

#include "numtracer/numeric/numeric_contract.hpp"
#include "numtracer/numeric/frame.hpp"
#include "engine_test_util.hpp"

using numtracer::Cx;
using namespace numtracer;
namespace nm = numtracer::numeric;
namespace nd = numtracer::numeric::ndetail;

// momenta: k = vid 0, a = vid 1, b = vid 2, q = vid 3 (a second projector momentum)
constexpr int kVid = 0, aVid = 1, bVid = 2, qVid = 3;
// atoms: 0 = 1/k², 1 = 1/q², 2 = 1/|k⃗|² (spatial, for E/M)
static LorentzFactor PT(int a, int b) { return ntest::fprojT(a, b, {{1.0, kVid}}, 0); }
static LorentzFactor PL(int a, int b) { return ntest::fprojL(a, b, {{1.0, kVid}}, 0); }
static LorentzFactor PE(int a, int b) { return ntest::fprojE(a, b, {{1.0, kVid}}, 0, 2); }
static LorentzFactor PM(int a, int b) { return ntest::fprojM(a, b, {{1.0, kVid}}, 2); }
static LorentzFactor PTq(int a, int b) { return ntest::fprojT(a, b, {{1.0, qVid}}, 1); }
static LorentzFactor vec(int id, int vid) { return ntest::fvec(id, {{1.0, vid}}); }

static int fails = 0;
static void check(bool ok, const char *name)
{
  std::printf("  %-44s %s\n", name, ok ? "ok" : "FAIL");
  if (!ok) ++fails;
}

// run fuse_projectors on a copy and report the folded list + coeff
struct Fold {
  std::vector<LorentzFactor> e;
  Cx c;
};
static Fold fold(std::vector<LorentzFactor> e, Cx c = {1, 0})
{
  nd::fuse_projectors(e, c);
  return {std::move(e), c};
}
// does the folded list hold exactly one projector of `kind` with the index set {x,y}?
static bool oneProj(const Fold &f, LorentzFactor::Kind kind, int x, int y)
{
  if (f.e.size() != 1 || f.e[0].kind != kind) return false;
  const int a = f.e[0].a, b = f.e[0].b;
  return (a == x && b == y) || (a == y && b == x);
}

int main()
{
  std::printf("== projector fusion: structural ==\n");

  // ---- idempotency, single shared inner index 12, ALL four symmetric slot variants -> P(10,11) ----
  check(oneProj(fold({PT(10, 12), PT(12, 11)}), LorentzFactor::ProjT, 10, 11), "PT(10,12).PT(12,11) -> PT(10,11)");
  check(oneProj(fold({PT(12, 10), PT(12, 11)}), LorentzFactor::ProjT, 10, 11), "PT(12,10).PT(12,11) -> PT(10,11)");
  check(oneProj(fold({PT(10, 12), PT(11, 12)}), LorentzFactor::ProjT, 10, 11), "PT(10,12).PT(11,12) -> PT(10,11)");
  check(oneProj(fold({PT(12, 10), PT(11, 12)}), LorentzFactor::ProjT, 10, 11), "PT(12,10).PT(11,12) -> PT(10,11)");

  // ---- idempotency for L / E / M ----
  check(oneProj(fold({PL(10, 12), PL(12, 11)}), LorentzFactor::ProjL, 10, 11), "PL.PL -> PL");
  check(oneProj(fold({PE(10, 12), PE(12, 11)}), LorentzFactor::ProjE, 10, 11), "PE.PE -> PE");
  check(oneProj(fold({PM(10, 12), PM(12, 11)}), LorentzFactor::ProjM, 10, 11), "PM.PM -> PM");

  // ---- a chain of three collapses to one (inner 12,13 are dummies; 10,11 survive) ----
  check(oneProj(fold({PT(10, 12), PT(12, 13), PT(13, 11)}), LorentzFactor::ProjT, 10, 11), "PT chain x3 -> PT(10,11)");

  // ---- full contraction (both indices shared) -> scalar trace; element list emptied ----
  {
    Fold f = fold({PT(10, 11), PT(10, 11)});
    check(f.e.empty() && f.c.re == 3 && f.c.im == 0, "PT.PT (both shared) -> tr=3");
  }
  {
    Fold f = fold({PT(10, 11), PT(11, 10)}); // swapped indices, still a full trace
    check(f.e.empty() && f.c.re == 3, "PT(10,11).PT(11,10) -> tr=3");
  }
  check(fold({PL(10, 11), PL(10, 11)}).c.re == 1, "PL.PL (both shared) -> tr=1");
  check(fold({PE(10, 11), PE(10, 11)}).c.re == 1, "PE.PE (both shared) -> tr=1");
  check(fold({PM(10, 11), PM(10, 11)}).c.re == 2, "PM.PM (both shared) -> tr=2");

  // ---- orthogonality, single shared index -> term vanishes (coeff 0) ----
  check(fold({PL(10, 12), PT(12, 11)}).c.re == 0, "PL.PT -> 0");
  check(fold({PE(10, 12), PM(12, 11)}).c.re == 0, "PE.PM -> 0");
  check(fold({PL(10, 12), PE(12, 11)}).c.re == 0, "PL.PE -> 0");
  check(fold({PL(10, 12), PM(12, 11)}).c.re == 0, "PL.PM -> 0");
  // orthogonality, both indices shared (a trace of orthogonal projectors) -> 0
  check(fold({PL(10, 11), PT(10, 11)}).c.re == 0, "PL.PT (both shared) -> 0");

  // ---- negatives: must NOT fold (list unchanged) ----
  auto unchanged = [](const Fold &f, std::size_t n) { return f.e.size() == n && f.c.re == 1; };
  check(unchanged(fold({PT(10, 12), PTq(12, 11)}), 2), "different momenta: not folded");
  check(unchanged(fold({PT(10, 12), PE(12, 11)}), 2), "T.E same momentum: left untouched");
  check(unchanged(fold({PT(10, 12), PM(12, 11)}), 2), "T.M same momentum: left untouched");
  // shared index 12 is NOT a clean dummy (also on a vector) -> no fold
  check(unchanged(fold({PT(10, 12), PT(12, 11), vec(12, aVid)}), 3), "non-dummy shared index: not folded");
  // single projector -> early-exit, unchanged
  check(unchanged(fold({PT(10, 11)}), 1), "single projector: no-op");

  // ───────────────────────── (B) value-preserving end-to-end via contract ─────────────────────────
  std::printf("== projector fusion: value-preserving ==\n");
  constexpr int nsym = 16;
  nm::Frame env(ntest::names(nsym));
  std::vector<std::array<nm::Poly, 4>> comp(4);
  for (int mu = 0; mu < 4; ++mu) {
    comp[kVid][(std::size_t)mu] = env.var(0 + mu);
    comp[aVid][(std::size_t)mu] = env.var(4 + mu);
    comp[bVid][(std::size_t)mu] = env.var(8 + mu);
    comp[qVid][(std::size_t)mu] = env.var(12 + mu);
  }
  nm::Poly k2 = env.zero(), ks2 = env.zero();
  for (int mu = 0; mu < 4; ++mu)
    k2 = k2 + comp[kVid][(std::size_t)mu] * comp[kVid][(std::size_t)mu];
  for (int mu = 1; mu < 4; ++mu)
    ks2 = ks2 + comp[kVid][(std::size_t)mu] * comp[kVid][(std::size_t)mu];
  const std::vector<nm::Poly> atomDen = {k2, env.zero(), ks2}; // atom 1 (q²) unused here

  auto NV = [&](std::initializer_list<LorentzFactor> e) {
    return ntest::contract(env, {}, LorentzNet{LorentzTerm{Cx{1, 0}, std::vector<LorentzFactor>(e)}}, comp, atomDen);
  };
  // redundant (fused internally) vs hand-simplified — must be equal MPolys
  const nm::Poly redT = NV({vec(10, aVid), PT(10, 12), PT(12, 11), vec(11, bVid)});
  const nm::Poly simT = NV({vec(10, aVid), PT(10, 11), vec(11, bVid)});
  const nm::Poly redChain = NV({vec(10, aVid), PT(10, 12), PT(12, 13), PT(13, 11), vec(11, bVid)});
  const nm::Poly orthoLT = NV({vec(10, aVid), PL(10, 12), PT(12, 11), vec(11, bVid)});
  const nm::Poly trClosed = NV({PT(10, 11), PT(10, 11)});

  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  double wEq = 0, wChain = 0, wOrtho = 0, wTr = 0;
  for (int it = 0; it < 200; ++it) {
    std::vector<double> x(nsym);
    for (double &v : x) v = U(rng);
    const double kk = nm::ndetail::eval(k2, x, {}).re;
    const double kss = nm::ndetail::eval(ks2, x, {}).re;
    const std::vector<double> av = {1.0 / kk, 0.0, 1.0 / kss};
    auto e = [&](const nm::Poly &p) { return nm::ndetail::eval(p, x, av).re; };
    wEq = std::max(wEq, std::fabs(e(redT) - e(simT)));
    wChain = std::max(wChain, std::fabs(e(redChain) - e(simT)));
    wOrtho = std::max(wOrtho, std::fabs(e(orthoLT)));
    wTr = std::max(wTr, std::fabs(e(trClosed) - 3.0));
  }
  check(wEq < 1e-12, "a.PT.PT.b == a.PT.b");
  check(wChain < 1e-12, "a.PT.PT.PT.b == a.PT.b");
  check(wOrtho < 1e-12, "a.PL.PT.b == 0");
  check(wTr < 1e-12, "PT.PT (closed) == tr PT = 3");

  std::printf(fails ? "\nTESTS FAILED (%d)\n" : "\nALL TESTS PASSED\n", fails);
  return fails ? 1 : 0;
}
