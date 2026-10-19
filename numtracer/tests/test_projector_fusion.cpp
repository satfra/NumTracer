/// @file test_projector_fusion.cpp
/// @brief Unit tests for the projector-algebra fusion (`ndetail::fuse_projectors`): two projectors on
///        the SAME momentum contracting over a shared dummy index collapse by the algebra —
///        idempotency `P·P→P` (one shared index) / `tr P` (both shared: T→3,L→1,E→1,M→2), and
///        orthogonality `L⟂{T,E,M}`, `E⟂M` → the term vanishes, and absorption `T·E→E`, `T·M→M`
///        (E and M live inside the transverse subspace, so P_T acts there as the identity).
///
/// Two layers:
///   (A) STRUCTURAL — call fuse_projectors on a hand-built element list and assert the folded list +
///       coeff (deterministic; pins the fold logic, incl. all symmetric index-slot variants, chains,
///       traces, orthogonality, and the no-op negatives).
///   (B) VALUE — contract a redundant net vs its hand-simplified form through the contraction and check
///       they agree at random points (fusion is value-preserving).
///   (C) `mul` drops a product term outright when it pairs two orthogonal projectors on one momentum
///       over a shared index — before the term is ever built in full. Structural checks here; the
///       value checks sit in (B).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <vector>

#include "numtracer/engine/contract.hpp"
#include "numtracer/engine/frame.hpp"
#include "engine_test_util.hpp"

using numtracer::Cx;
using namespace numtracer;
namespace nm = numtracer;
namespace nd = numtracer::engine::ndetail;

// momenta: k = vid 0, a = vid 1, b = vid 2, q = vid 3 (a second projector momentum)
constexpr int kVid = 0, aVid = 1, bVid = 2, qVid = 3;
// atoms: 0 = 1/k², 1 = 1/q², 2 = 1/|k⃗|² (spatial, for E/M)
static LorentzFactor PT(int a, int b) { return ntest::fprojT(a, b, {{1.0, kVid}}, 0); }
static LorentzFactor PL(int a, int b) { return ntest::fprojL(a, b, {{1.0, kVid}}, 0); }
static LorentzFactor PE(int a, int b) { return ntest::fprojE(a, b, {{1.0, kVid}}, 0, 2); }
static LorentzFactor PM(int a, int b) { return ntest::fprojM(a, b, {{1.0, kVid}}, 2); }
static LorentzFactor PTq(int a, int b) { return ntest::fprojT(a, b, {{1.0, qVid}}, 1); }
static LorentzFactor PMq(int a, int b) { return ntest::fprojM(a, b, {{1.0, qVid}}, 2); }
static LorentzFactor vec(int id, int vid) { return ntest::fvec(id, {{1.0, vid}}); }
static LorentzFactor met(int a, int b) { return ntest::fmet(a, b); }

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

  // ---- absorption T·E -> E, T·M -> M (either order); the result keeps the E/M factor's atoms ----
  auto keepsAtoms = [](const Fold &f, int atom, int atomS) { return f.e[0].atom == atom && f.e[0].atomS == atomS; };
  {
    Fold f = fold({PT(10, 12), PE(12, 11)});
    check(oneProj(f, LorentzFactor::ProjE, 10, 11) && keepsAtoms(f, 0, 2), "PT.PE -> PE(10,11)");
  }
  {
    Fold f = fold({PE(10, 12), PT(12, 11)});
    check(oneProj(f, LorentzFactor::ProjE, 10, 11) && keepsAtoms(f, 0, 2), "PE.PT -> PE(10,11)");
  }
  {
    Fold f = fold({PT(10, 12), PM(12, 11)});
    check(oneProj(f, LorentzFactor::ProjM, 10, 11) && keepsAtoms(f, -1, 2), "PT.PM -> PM(10,11)");
  }
  {
    Fold f = fold({PM(10, 12), PT(12, 11)});
    check(oneProj(f, LorentzFactor::ProjM, 10, 11) && keepsAtoms(f, -1, 2), "PM.PT -> PM(10,11)");
  }
  check(oneProj(fold({PT(10, 12), PE(12, 13), PT(13, 11)}), LorentzFactor::ProjE, 10, 11), "PT.PE.PT -> PE(10,11)");
  check(fold({PE(10, 12), PT(12, 13), PM(13, 11)}).c.re == 0, "PE.PT.PM -> 0");
  {
    Fold f = fold({PT(10, 11), PE(10, 11)});
    check(f.e.empty() && f.c.re == 1, "PT.PE (both shared) -> tr=1");
  }
  {
    Fold f = fold({PM(10, 11), PT(11, 10)});
    check(f.e.empty() && f.c.re == 2, "PM.PT (both shared) -> tr=2");
  }

  // ---- negatives: must NOT fold (list unchanged) ----
  auto unchanged = [](const Fold &f, std::size_t n) { return f.e.size() == n && f.c.re == 1; };
  check(unchanged(fold({PT(10, 12), PTq(12, 11)}), 2), "different momenta: not folded");
  check(unchanged(fold({PTq(10, 12), PE(12, 11)}), 2), "T.E different momenta: not folded");
  // shared index 12 is NOT a clean dummy (also on a vector) -> no fold
  check(unchanged(fold({PT(10, 12), PT(12, 11), vec(12, aVid)}), 3), "non-dummy shared index: not folded");
  // single projector -> early-exit, unchanged
  check(unchanged(fold({PT(10, 11)}), 1), "single projector: no-op");

  // ───────────────────────── (C) mul drops orthogonal projector pairs ─────────────────────────
  std::printf("== mul: orthogonal projector pairs ==\n");
  auto L1 = [](LorentzFactor f) { return leaf(std::move(f)); };
  auto zeroNet = [](const LorentzNet &n) { return n.size() == 1 && n[0].coeff.re == 0 && n[0].coeff.im == 0; };
  auto kindsOf = [](const LorentzTerm &t) {
    std::vector<int> k;
    for (const auto &f : t.e) k.push_back(f.kind);
    return k;
  };
  check(zeroNet(mul(L1(PE(10, 12)), L1(PM(12, 11)))), "E(10,12).M(12,11) -> explicit zero term");
  check(zeroNet(mul(L1(PM(10, 11)), L1(PE(11, 10)))), "M.E (both shared) -> explicit zero term");
  check(zeroNet(mul(L1(PL(10, 12)), L1(PT(12, 11)))), "L.T -> explicit zero term");
  check(zeroNet(mul(L1(vec(10, aVid)), L1(PE(10, 12)), L1(PM(12, 11)), L1(vec(11, bVid)))),
        "a.E.M.b -> stays zero through later muls");
  {
    const LorentzNet n = mul(add(L1(PE(10, 12)), L1(PM(10, 12))), add(L1(PE(12, 11)), L1(PM(12, 11))));
    const bool ok = n.size() == 2 && kindsOf(n[0]) == std::vector<int>{LorentzFactor::ProjE, LorentzFactor::ProjE} &&
                    kindsOf(n[1]) == std::vector<int>{LorentzFactor::ProjM, LorentzFactor::ProjM};
    check(ok, "(E+M).(E+M) -> E.E + M.M only");
  }
  check(mul(L1(PE(10, 12)), L1(PMq(12, 11))).size() == 1 && mul(L1(PE(10, 12)), L1(PMq(12, 11)))[0].e.size() == 2,
        "E(k).M(q): different momenta, kept");
  check(mul(L1(PE(10, 11)), L1(PM(12, 13)))[0].e.size() == 2, "E.M without a shared index, kept");
  check(mul(L1(PT(10, 12)), L1(PE(12, 11)))[0].e.size() == 2, "T.E (not orthogonal), kept");

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
  // Absorption, graded against references whose metric δ(12,13) blocks fusion, so they are
  // contracted from the dense tensors and do not depend on the rule under test.
  const nm::Poly redTE = NV({vec(10, aVid), PT(10, 12), PE(12, 11), vec(11, bVid)});
  const nm::Poly refTE = NV({vec(10, aVid), PT(10, 12), met(12, 13), PE(13, 11), vec(11, bVid)});
  const nm::Poly redMT = NV({vec(10, aVid), PM(10, 12), PT(12, 11), vec(11, bVid)});
  const nm::Poly refMT = NV({vec(10, aVid), PM(10, 12), met(12, 13), PT(13, 11), vec(11, bVid)});
  const nm::Poly trTE = NV({PT(10, 11), PE(10, 11)});
  const nm::Poly trTM = NV({PT(10, 11), PM(11, 10)});
  // (C) mul's pruning: a product that is entirely orthogonal contracts to 0 (not to the empty net's 1),
  // and a pruned sum matches the same sum written out term by term.
  auto CN = [&](const LorentzNet &n) { return ntest::contract(env, {}, n, comp, atomDen); };
  const nm::Poly zeroEM = CN(mul(L1(vec(10, aVid)), L1(PE(10, 12)), L1(PM(12, 11)), L1(vec(11, bVid))));
  const nm::Poly prunedEM =
      CN(mul(L1(vec(10, aVid)), add(L1(PE(10, 12)), L1(PM(10, 12))), add(L1(PE(12, 11)), L1(PM(12, 11))), L1(vec(11, bVid))));
  auto T4 = [&](LorentzFactor x, LorentzFactor y) {
    return LorentzTerm{Cx{1, 0}, {vec(10, aVid), std::move(x), std::move(y), vec(11, bVid)}};
  };
  const nm::Poly fullEM = CN(LorentzNet{T4(PE(10, 12), PE(12, 11)), T4(PE(10, 12), PM(12, 11)),
                                        T4(PM(10, 12), PE(12, 11)), T4(PM(10, 12), PM(12, 11))});

  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  double wEq = 0, wChain = 0, wOrtho = 0, wTr = 0, wTE = 0, wMT = 0, wTrTE = 0, wTrTM = 0, scale = 0;
  double wZeroEM = 0, wPruned = 0, scaleEM = 0;
  for (int it = 0; it < 200; ++it) {
    std::vector<double> x(nsym);
    for (double &v : x) v = U(rng);
    const double kk = nm::pdetail::eval(k2, x, {}).re;
    const double kss = nm::pdetail::eval(ks2, x, {}).re;
    const std::vector<double> av = {1.0 / kk, 0.0, 1.0 / kss};
    auto e = [&](const nm::Poly &p) { return nm::pdetail::eval(p, x, av).re; };
    wEq = std::max(wEq, std::fabs(e(redT) - e(simT)));
    wChain = std::max(wChain, std::fabs(e(redChain) - e(simT)));
    wOrtho = std::max(wOrtho, std::fabs(e(orthoLT)));
    wTr = std::max(wTr, std::fabs(e(trClosed) - 3.0));
    wTE = std::max(wTE, std::fabs(e(redTE) - e(refTE)));
    wMT = std::max(wMT, std::fabs(e(redMT) - e(refMT)));
    wTrTE = std::max(wTrTE, std::fabs(e(trTE) - 1.0));
    wTrTM = std::max(wTrTM, std::fabs(e(trTM) - 2.0));
    scale = std::max({scale, std::fabs(e(refTE)), std::fabs(e(refMT))});
    wZeroEM = std::max(wZeroEM, std::fabs(e(zeroEM)));
    wPruned = std::max(wPruned, std::fabs(e(prunedEM) - e(fullEM)));
    scaleEM = std::max(scaleEM, std::fabs(e(fullEM)));
  }
  check(wEq < 1e-12, "a.PT.PT.b == a.PT.b");
  check(wChain < 1e-12, "a.PT.PT.PT.b == a.PT.b");
  check(wOrtho < 1e-12, "a.PL.PT.b == 0");
  check(wTr < 1e-12, "PT.PT (closed) == tr PT = 3");
  // relative to the largest reference value: 1/|k⃗|² makes the samples range widely
  check(wTE < 1e-10 * std::max(1.0, scale), "a.PT.PE.b == a.PT.δ.PE.b");
  check(wMT < 1e-10 * std::max(1.0, scale), "a.PM.PT.b == a.PM.δ.PT.b");
  check(wTrTE < 1e-10, "PT.PE (closed) == tr PE = 1");
  check(wTrTM < 1e-10, "PT.PM (closed) == tr PM = 2");
  check(wZeroEM == 0, "a.E.M.b contracts to exactly 0");
  check(wPruned < 1e-10 * std::max(1.0, scaleEM), "pruned a.(E+M).(E+M).b == all four terms");

  // A Dirac loop whose free indices were closed only by the pruned projectors: the explicit zero term
  // that mul leaves behind must contract to 0, not trip the open-index guard on the loop's legs.
  {
    bool threw = false, zero = false;
    try {
      const nm::Poly p = ntest::contract(env, DiracChain{gamma(LorentzIndex{10}), gamma(LorentzIndex{11})},
                                         mul(L1(PE(10, 12)), L1(PM(12, 11))), comp, atomDen);
      zero = p.terms.empty();
    } catch (const std::exception &) {
      threw = true;
    }
    check(!threw && zero, "tr(γ^μ γ^ν).E.M (all pruned) == 0, no throw");
  }

  std::printf(fails ? "\nTESTS FAILED (%d)\n" : "\nALL TESTS PASSED\n", fails);
  return fails ? 1 : 0;
}
