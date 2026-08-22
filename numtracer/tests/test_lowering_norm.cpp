/// @file test_lowering_norm.cpp
/// @brief Differential test for the lowering strategies in `codegen/lower.hpp`.
///
/// The lowering is the one stage `ctest` is otherwise blind to without invoking Wolfram: every other
/// C++ test compiles a PRE-GENERATED kernel. This test drives @ref numtracer::network::horner
/// directly on hand-built monomial sets, lowers each one BOTH ways (scalar-normalised and plain),
/// interprets each instruction stream, and checks it against a direct evaluation of the polynomial.
/// It needs no Wolfram, no CUDA, and runs in well under a second.
///
/// What it is really guarding is the aliasing of `-1`: to the SSA builder a negative slot is
/// structural ZERO, while to an @ref numtracer::network::NVal shape it is the structural ONE. Getting
/// that wrong does not crash — it silently emits a kernel that returns `1.0`, or `0.0`, instead of the
/// polynomial. The `pure constant`, `cancels to zero` and `empty` cases below are exactly those.

#include "numtracer/codegen/gen.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace numtracer::network;
using namespace numtracer::network::rdetail;

namespace
{
  constexpr int kNumVars = 10;

  const bool kModes[] = {false, true};
  const char *kModeNames[] = {"plain", "normalised"};
  constexpr int kNumModes = 2;

  /// Evaluate an emitted instruction stream — the semantics `emit_stmt` prints.
  double interpret(const std::vector<RInstr> &ins, int root, const double *f)
  {
    if (root < 0) return 0.0;
    std::vector<double> v(ins.size(), 0.0);
    auto val = [&](int r) { return r < 0 ? 0.0 : v[static_cast<std::size_t>(r)]; };
    for (std::size_t i = 0; i < ins.size(); ++i) {
      const RInstr &in = ins[i];
      switch (in.op) {
      case RCONST: v[i] = in.value; break;
      case RVAR: v[i] = f[in.a]; break;
      case RADD: v[i] = val(in.a) + val(in.b); break;
      case RMUL: v[i] = val(in.a) * val(in.b); break;
      default: v[i] = -val(in.a); break;
      }
    }
    return v[static_cast<std::size_t>(root)];
  }

  /// Evaluate the monomial set directly — the reference every lowering must reproduce.
  double evaluate(const std::vector<LMono> &ms, const double *f)
  {
    double acc = 0;
    for (const LMono &m : ms) {
      double t = m.c;
      for (auto [id, e] : m.vp)
        for (int k = 0; k < e; ++k) t *= f[id];
      acc += t;
    }
    return acc;
  }

  int failures = 0;

  /// Lower @p ms under every mode and check each against @ref evaluate at @p f.
  void check(const std::string &name, const std::vector<LMono> &ms, const double *f, double tol = 1e-11)
  {
    const double ref = evaluate(ms, f);
    for (int m = 0; m < kNumModes; ++m) {
      RBuilder b;
      const int root = scale_into(b, horner(b, ms, kModes[m]));
      const double got = interpret(b.ins, root, f);
      const double scale = std::max({std::fabs(ref), std::fabs(got), 1e-300});
      if (!(std::fabs(got - ref) / scale < tol)) {
        std::printf("  FAIL %-26s [%s] ref=%.17g got=%.17g\n", name.c_str(), kModeNames[m], ref, got);
        ++failures;
      }
    }
  }

  LMono mono(double c, std::vector<std::pair<int, int>> vp)
  {
    LMono m;
    m.c = c;
    m.vp = std::move(vp);
    std::sort(m.vp.begin(), m.vp.end());
    return m;
  }
} // namespace

int main()
{
  double f[kNumVars];
  std::mt19937 rng(20260821);
  std::uniform_real_distribution<double> ud(0.3, 2.5);
  for (int i = 0; i < kNumVars; ++i) f[i] = ud(rng);

  // ---- the corner cases where the zero/one aliasing bites -------------------------------------
  check("empty", {}, f);
  check("pure constant", {mono(-7.25, {})}, f);
  check("constant zero", {mono(0.0, {})}, f);
  check("cancels to zero", {mono(3.5, {{0, 1}}), mono(-3.5, {{0, 1}})}, f);
  check("coefficient +1", {mono(1.0, {{0, 1}, {1, 2}})}, f);
  check("coefficient -1", {mono(-1.0, {{0, 1}, {1, 2}})}, f);
  // exact +-1 ratios: `rmul` folds x(-1) to an RNEG, which the printer then spells as a subtraction
  check("ratio +1", {mono(4.0, {{0, 1}}), mono(4.0, {{1, 1}})}, f);
  check("ratio -1", {mono(4.0, {{0, 1}}), mono(-4.0, {{1, 1}})}, f);
  // the relative noise prune bounds the coefficient spread at ~1e9; go to the edge of it
  check("1e9 coefficient spread", {mono(1e9, {{0, 1}}), mono(1.0, {{1, 1}}), mono(-1e-9, {{2, 1}})}, f, 1e-10);
  check("high powers", {mono(2.0, {{0, 5}}), mono(-3.0, {{0, 3}, {1, 2}}), mono(1.5, {{1, 4}})}, f);

  // ---- the shape scalar normalisation exists for: few shapes, many distinct scalars -------------
  // (a plain Horner pushes the scalar into the leaves, so each of these interns separately)
  {
    const std::vector<std::pair<int, int>> shapes[] = {
        {{0, 1}, {1, 1}}, {{2, 1}}, {{0, 1}, {2, 2}}, {{1, 1}, {3, 1}}};
    const double scalars[] = {295.60333782509002, -886.81001347527001, 443.40500673763, 36.950417228135997,
                              -32.0,              1.0,                 -1.0,            110.85125168441};
    std::vector<LMono> ms;
    for (const auto &sh : shapes)
      for (double c : scalars) ms.push_back(mono(c, sh));
    check("shared shapes x scalars", ms, f);
  }

  // ---- a shared COMPOSITE factor: the pivot only ever factors out an env VARIABLE, so this stays
  // distributed. Kept as a shape the lowering must still get RIGHT, not as one it optimises. -------
  check("shared composite factor",
        {mono(2.0, {{0, 1}, {1, 1}, {2, 1}}), mono(3.0, {{0, 1}, {1, 1}, {3, 1}, {4, 1}}),
         mono(-5.0, {{0, 1}, {1, 1}, {5, 2}})},
        f);

  // ---- randomised sweep ------------------------------------------------------------------------
  long slots[kNumModes] = {0, 0};
  for (int trial = 0; trial < 500; ++trial) {
    const int nshape = 1 + static_cast<int>(rng() % 6);
    std::vector<std::vector<std::pair<int, int>>> shapes;
    for (int i = 0; i < nshape; ++i) {
      std::vector<std::pair<int, int>> vp;
      const int nf = static_cast<int>(rng() % 4);
      for (int j = 0; j < nf; ++j) vp.push_back({static_cast<int>(rng() % kNumVars), 1 + static_cast<int>(rng() % 3)});
      std::sort(vp.begin(), vp.end());
      vp.erase(std::unique(vp.begin(), vp.end(), [](auto &a, auto &b) { return a.first == b.first; }), vp.end());
      shapes.push_back(std::move(vp));
    }
    const double base[] = {1.0, -1.0, 2.0, 295.60333782509002, -886.81001347527001, 0.5, 1e-7, 3.0};
    std::vector<LMono> ms;
    const int nm = 1 + static_cast<int>(rng() % 40);
    for (int i = 0; i < nm; ++i)
      ms.push_back(mono(base[rng() % 8] * static_cast<double>(1 + static_cast<int>(rng() % 4)),
                        shapes[rng() % shapes.size()]));
    check("random " + std::to_string(trial), ms, f);

    for (int m = 0; m < kNumModes; ++m) {
      RBuilder b;
      scale_into(b, horner(b, ms, kModes[m]));
      slots[m] += static_cast<long>(b.ins.size());
    }
  }

  std::printf("lowering A/B: %d failure(s); random-sweep slots plain=%ld normalised=%ld\n",
              failures, slots[0], slots[1]);
  if (failures == 0) std::printf("TESTS PASSED\n");
  return failures == 0 ? 0 : 1;
}
