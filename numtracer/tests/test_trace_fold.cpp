// The generator's sub-term dedup (numeric/trace_fold.hpp): phase A contracts each DISTINCT Dirac
// trace once, phase B folds each net's traces with its own scalars as a balanced binary-counter tree.
//
// Why this test exists. The dedup only does anything on a flow with MANY sub-terms per net, and every
// committed flow in tests/gen has exactly ONE — so the tree fold degenerates to a single push there
// and the whole suite, codegen test included, is blind to it. The fold is also the piece most likely
// to be silently wrong: a carry-propagation stack that drops or double-counts a partial still returns
// a plausible polynomial. So exercise it directly, against the naive left fold it replaces, at the
// sizes and cache regimes the real flows hit:
//
//   - fold_net == the left fold  (reassociation must not change the VALUE; MPoly/DPoly addition is
//     exact term collection over Cx, so this is an exact equality, not a tolerance)
//   - every cache regime: nCache = 0 (nothing resident, all recomputed), 1, n/2, n (all resident)
//   - a trace beyond nCache is recomputed EXACTLY ONCE (the singleton-eviction assumption: Codegen.m
//     only leaves refCount==1 traces uncached, so a recompute must never be duplicated work)
//   - both backends: MPoly (plain) and DPoly (dressed)
//   - the degenerate shapes that the binary counter gets wrong if the carry is off by one: 0, 1, 2, 3
//     and 7 terms (7 = 0b111 leaves three unmerged partials on the stack -> the leftover sweep)
#include "numtracer/numeric/env.hpp"
#include "numtracer/numeric/trace_fold.hpp"

#include <cstdio>
#include <numeric>
#include <vector>

using namespace numtracer;
using namespace numtracer::numeric;

namespace
{
  int failures = 0;

  void check(bool ok, const char *what)
  {
    if (!ok) {
      std::printf("  FAIL: %s\n", what);
      ++failures;
    }
  }

  constexpr int NSYM = 3;

  // A distinct, non-trivial polynomial per trace id: c(k)·x0^(k%3) · x1^(k%2) + k·x2. Distinct enough
  // that a dropped or double-counted term cannot cancel out of the total by luck.
  MPoly trace_of(int k)
  {
    LorentzEnv env(NSYM);
    MPoly p = env.mono({k % 3, k % 2, 0}, Cx{1.0 + k, 0.5 * k});
    return p + env.mono({0, 0, 1}, Cx{static_cast<double>(k), 0.0});
  }

  bool same(const MPoly &a, const MPoly &b)
  {
    if (a.terms.size() != b.terms.size()) return false;
    for (std::size_t i = 0; i < a.terms.size(); ++i) {
      if (!(a.terms[i].first == b.terms[i].first)) return false;
      if (a.terms[i].second.re != b.terms[i].second.re || a.terms[i].second.im != b.terms[i].second.im) return false;
    }
    return true;
  }

  bool same(const DPoly &a, const DPoly &b)
  {
    if (a.terms.size() != b.terms.size()) return false;
    for (std::size_t i = 0; i < a.terms.size(); ++i) {
      if (a.terms[i].first != b.terms[i].first) return false;
      if (!same(a.terms[i].second, b.terms[i].second)) return false;
    }
    return true;
  }

  // A stand-in polynomial that records how many ELEMENTS each addition touches. `size` models the
  // term count; the cost of adding two polynomials is linear in their combined size (MPoly::operator+
  // is a merge). Summing that over the fold is the quantity the tree fold is built to keep at
  // O(n log n) — a left fold against a growing accumulator makes it O(n²).
  struct Counted {
    long size = 0;
    static long touches;
    Counted() = default;
    explicit Counted(int) {}                                  // the P(nsym) empty accumulator
    static Counted terms(long s) { Counted c; c.size = s; return c; }
  };
  long Counted::touches = 0;

  Counted operator+(const Counted &a, const Counted &b)
  {
    Counted::touches += a.size + b.size;                      // a merge touches both operands
    return Counted::terms(a.size + b.size);
  }
  Counted scale_trace(int, const Counted &p, Cx) { return p; } // found by ADL from inside fold_net

  // What fold_net replaces: Σ_j sc[j]·trace(idx[j]), left-folded against one growing accumulator.
  template <class P, class TraceFn>
  P left_fold(const std::vector<int> &idx, const std::vector<Cx> &sc, TraceFn &&trace)
  {
    P acc = zero_like<P>(NSYM); // MPoly/DPoly route through the private-factory attorney; Counted uses its ctor
    for (std::size_t j = 0; j < idx.size(); ++j)
      acc = acc + scale_trace(NSYM, trace(idx[j]), sc[j]);
    return acc;
  }

  // One (backend, net shape, cache regime) case.
  template <class P, class MakeTrace>
  void run_case(const char *tag, std::size_t n, long nCache, MakeTrace &&mk)
  {
    std::vector<int> idx(n);
    std::vector<Cx> sc(n);
    for (std::size_t j = 0; j < n; ++j) {
      idx[j] = static_cast<int>(j);            // every trace distinct => refCount 1 => recompute is 1x
      sc[j] = Cx{1.0 + 0.25 * j, -0.125 * j};  // non-uniform scalars: a swapped scalar shows up
    }

    // Phase A over the resident prefix, exactly as the generated main does it.
    std::vector<P> T = contract_traces<P>(NSYM, nCache, /*W=*/4, mk);
    check(T.size() == static_cast<std::size_t>(nCache < 0 ? 0 : nCache), "phase A table size");

    // Count recomputes so we can prove nothing beyond the cache is computed twice.
    std::vector<int> hits(n, 0);
    auto counting = [&](int k) {
      if (static_cast<std::size_t>(k) < hits.size()) ++hits[static_cast<std::size_t>(k)];
      return mk(k);
    };

    const P got = fold_net<P>(NSYM, idx, sc, T, nCache, counting);
    const P want = left_fold<P>(idx, sc, mk);

    char msg[160];
    std::snprintf(msg, sizeof msg, "%s n=%zu nCache=%ld: tree fold == left fold", tag, n, nCache);
    check(same(got, want), msg);

    for (std::size_t k = 0; k < n; ++k) {
      const int expect = (static_cast<long>(k) < nCache) ? 0 : 1; // resident => 0 recomputes, else 1
      if (hits[k] != expect) {
        std::snprintf(msg, sizeof msg, "%s n=%zu nCache=%ld: trace %zu recomputed %dx (want %dx)", tag, n,
                      nCache, k, hits[k], expect);
        check(false, msg);
        break;
      }
    }
  }

  template <class P, class MakeTrace> void run_backend(const char *tag, MakeTrace &&mk)
  {
    // 0/1/2/3 are the binary-counter edge cases; 7 = 0b111 leaves three partials for the leftover
    // sweep; 1000 is a heavy net, where the left fold's O(n^2) is what the tree is there to avoid.
    for (std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{2}, std::size_t{3}, std::size_t{7},
                          std::size_t{64}, std::size_t{1000}}) {
      for (long nCache : {0L, 1L, static_cast<long>(n) / 2, static_cast<long>(n)}) {
        if (nCache > static_cast<long>(n)) continue;
        run_case<P>(tag, n, nCache, mk);
      }
    }
  }
} // namespace

int main()
{
  std::printf("trace_fold: tree fold vs left fold, all cache regimes\n");

  run_backend<MPoly>("MPoly", trace_of);
  run_backend<DPoly>("DPoly", [](int k) {
    // A dressed trace: the same kinematic polynomial under a dressing monomial that varies with k, so
    // the DPoly merge (which is keyed on the dressing monomial) is exercised rather than bypassed.
    LorentzEnv env(NSYM);
    DPoly d = env.dzero();
    d.add(dmono_sorted(DMono{k % 2, k % 3}), trace_of(k));
    return d;
  });

  // A net that references the SAME trace repeatedly — the shape the dedup actually produces (a hot
  // trace is referenced 2-32x across the flow). Left fold and tree fold must still agree.
  {
    const std::vector<int> idx = {2, 2, 5, 2, 5, 0, 5, 5};
    const std::vector<Cx> sc = {{1, 0}, {-1, 0.5}, {2, 0}, {0.5, -0.25}, {1, 1}, {3, 0}, {-2, 0}, {0.25, 0}};
    std::vector<MPoly> T = contract_traces<MPoly>(NSYM, 6, 4, trace_of);
    check(same(fold_net<MPoly>(NSYM, idx, sc, T, 6, trace_of), left_fold<MPoly>(idx, sc, trace_of)),
          "repeated-trace net: tree fold == left fold");
  }

  // The tree must actually be a TREE. Every correct fold adds each term exactly once, so a botched
  // carry that merges only one level still returns the RIGHT VALUE — it just degenerates back into
  // the left fold this exists to replace. A value check cannot see that (verified: turning the carry
  // loop into a single `if` passes every check above). So bound the WORK instead: count the operand
  // sizes touched across all additions, which is what a left fold makes quadratic.
  {
    constexpr std::size_t n = 1024;
    std::vector<int> idx(n);
    std::vector<Cx> sc(n, Cx{1, 0});
    for (std::size_t j = 0; j < n; ++j)
      idx[j] = static_cast<int>(j);

    const std::vector<Counted> none;
    Counted::touches = 0;
    (void)fold_net<Counted>(NSYM, idx, sc, none, 0, [](int) { return Counted::terms(1); });
    const long tree = Counted::touches;

    Counted::touches = 0;
    (void)left_fold<Counted>(idx, sc, [](int) { return Counted::terms(1); });
    const long left = Counted::touches;

    // n·log2(n) ≈ 10k vs n²/2 ≈ 524k. The 4x margin keeps this from being brittle while still being
    // orders of magnitude away from a left fold.
    char msg[160];
    std::snprintf(msg, sizeof msg, "tree fold is O(n log n): %ld touches vs left fold %ld (n=%zu)", tree,
                  left, n);
    check(tree * 4 < left, msg);
  }

  // STREAMING PHASE B (fold_groups_streaming) must be bit-identical to what it replaces: fold_nets
  // over every net, then an eager left fold of each group's members. It exists purely to bound RAM —
  // holding `window` group accumulators instead of every net polynomial — so any change to the VALUE
  // would be a bug, not a tradeoff. Three things are checked, because the first two can pass while the
  // third silently rots:
  //   - the accumulated polynomial per group is exactly equal (not within a tolerance)
  //   - the sink sees gi = 0,1,2,... strictly ascending, exactly once each. This is what keeps
  //     GlobalEnv intern order and the shared CSE instruction stream identical downstream, and it is
  //     NOT implied by the values being right — a batch drained out of order still sums correctly.
  //   - it holds across window = 1 (the RAM floor), an interior window that does not divide the group
  //     count evenly, and window = ngroups (which reproduces the old schedule exactly).
  {
    constexpr std::size_t nNet = 37;
    constexpr int NPT = 11; // distinct traces
    std::vector<std::vector<int>> sidx(nNet);
    std::vector<std::vector<Cx>> sc(nNet);
    for (std::size_t d = 0; d < nNet; ++d) {
      const std::size_t nsub = 1 + (d % 5); // skewed nets, like the real flows
      for (std::size_t j = 0; j < nsub; ++j) {
        sidx[d].push_back(static_cast<int>((d * 7 + j * 3) % NPT));
        sc[d].push_back(Cx{1.0 + 0.5 * static_cast<double>(d), -0.25 * static_cast<double>(j)});
      }
    }

    // A partition of [0, nNet) into ragged groups — the shape Codegen.m's GatherBy produces.
    std::vector<std::vector<int>> groups;
    for (std::size_t d = 0; d < nNet;) {
      const std::size_t k = 1 + (groups.size() % 4);
      std::vector<int> grp;
      for (std::size_t j = 0; j < k && d < nNet; ++j, ++d)
        grp.push_back(static_cast<int>(d));
      groups.push_back(std::move(grp));
    }

    // Per-net colour weights, and the scale the emitted DPoly branch uses.
    std::vector<Cx> colv(nNet);
    for (std::size_t d = 0; d < nNet; ++d)
      colv[d] = Cx{0.75 - 0.03 * static_cast<double>(d), 0.4 * static_cast<double>(d % 3)};

    for (long nCache : {0L, 5L, static_cast<long>(NPT)}) {
      const std::vector<MPoly> T = contract_traces<MPoly>(NSYM, nCache, /*W=*/4, trace_of);

      // Reference: exactly the code path being replaced.
      const std::vector<MPoly> mp = fold_nets<MPoly>(NSYM, sidx, sc, T, nCache, /*W=*/4, trace_of);
      std::vector<MPoly> want;
      for (const auto &grp : groups) {
        LorentzEnv env(NSYM);
        MPoly acc = env.zero();
        for (int d : grp)
          acc = acc + mp[static_cast<std::size_t>(d)] * env.constant(colv[static_cast<std::size_t>(d)]);
        want.push_back(std::move(acc));
      }

      for (long window : {1L, 5L, static_cast<long>(groups.size())}) {
        std::vector<MPoly> got;
        std::vector<std::size_t> order;
        LorentzEnv env(NSYM);
        env.fold_groups_streaming<MPoly>(
            sidx, sc, groups, T, nCache, /*W=*/4, window, trace_of,
            [&](int d, MPoly &&m) { return m * env.constant(colv[static_cast<std::size_t>(d)]); },
            [&](std::size_t gi, MPoly &&a) {
              order.push_back(gi);
              got.push_back(std::move(a));
            });

        char msg[192];
        std::snprintf(msg, sizeof msg, "streaming nCache=%ld window=%ld: one sink call per group",
                      nCache, window);
        check(got.size() == groups.size() && order.size() == groups.size(), msg);

        bool ascending = true;
        for (std::size_t i = 0; i < order.size(); ++i)
          if (order[i] != i) ascending = false;
        std::snprintf(msg, sizeof msg, "streaming nCache=%ld window=%ld: sink order strictly ascending",
                      nCache, window);
        check(ascending, msg);

        bool eq = got.size() == want.size();
        for (std::size_t i = 0; eq && i < got.size(); ++i)
          if (!same(got[i], want[i])) eq = false;
        std::snprintf(msg, sizeof msg, "streaming nCache=%ld window=%ld: == fold_nets + eager group fold",
                      nCache, window);
        check(eq, msg);
      }
    }

    // The partition checker must actually fire — it is the guard against a future Codegen.m grouping
    // change silently dropping or double-folding a net, and a checker that never reports is no guard.
    check_group_partition(groups, static_cast<long>(nNet)); // clean: does not exit

    // ...and it must fire on each way a grouping can be wrong. Assert on the PREDICATE, not on the
    // process: check_group_partition is fatal by design (a non-partition is a silently wrong kernel),
    // so it cannot be called with bad input from inside a test.
    {
      long d = 0, g = 0, o = 0;
      check(group_partition_stats({{0, 1}, {2}}, 3, d, g, o) && d == 0 && g == 0 && o == 0,
            "partition: exact cover is clean");
      check(!group_partition_stats({{0, 1}, {1, 2}}, 3, d, g, o) && d == 1,
            "partition: a net folded TWICE is caught");
      check(!group_partition_stats({{0}, {2}}, 3, d, g, o) && g == 1,
            "partition: a net DROPPED from the kernel is caught");
      check(!group_partition_stats({{0}, {5}}, 3, d, g, o) && o == 1,
            "partition: an out-of-range net index is caught");
      // A gap and a duplicate at once must report both, not stop at the first.
      check(!group_partition_stats({{0, 0}}, 2, d, g, o) && d == 1 && g == 1,
            "partition: duplicate AND gap both reported");
    }
  }

  // poly_bytes must actually track size — it is what the generator reports as the trace table's RAM
  // cost, and a constant would silently hide the one way this design can regress (peak RSS).
  LorentzEnv env(NSYM);
  check(poly_bytes(env.zero()) == 0, "poly_bytes: empty is 0");
  check(poly_bytes(trace_of(4)) > poly_bytes(env.constant(Cx{1, 0})), "poly_bytes: grows with terms");

  if (failures) {
    std::printf("trace_fold: %d FAILURES\n", failures);
    return 1;
  }
  std::printf("trace_fold: all checks passed\n");
  return 0;
}
