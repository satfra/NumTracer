/// @file trace_fold.hpp
/// @brief The generator's two contraction phases: contract each DISTINCT Dirac trace once (phase A),
///        then fold each net's traces with its own scalars (phase B).
///
/// The generator's cost is one trace contraction per `(net, sub-term)`, but the same
/// `(dnet, lnet, dch, dsl)` tuple recurs across nets and colour branches, so most of those
/// contractions recompute a trace that was already computed. Measured on dense flows: 30,807
/// contractions for 6,041 distinct traces (5.1x), and 246,456 for 32,784 (7.5x). Codegen.m therefore
/// emits a table of the distinct traces and, per net, the indices into it; this header contracts and
/// folds that table. Two independent wins:
///
///  - **the contraction shrinks by the redundancy factor** — each distinct trace is contracted once;
///  - **the parallel phase becomes a FLAT list of uniform work items.** Scheduling per *net* could
///    not use the machine: sub-terms per net are wildly skewed (max 2880 vs a median of 27), so the
///    single biggest net alone exceeded the ideal per-thread load and pinned utilisation at ~33%
///    however many cores were available.
///
/// Phase B reduces each net as a **balanced binary-counter tree** rather than a left fold against a
/// growing accumulator: heavy nets fold >1000 terms, and a left fold pays O(|acc|) on every one of
/// them. Reassociating perturbs the emitted kernel only in the last ulp (measured: worst relative
/// deviation ~7e-15), which is within the accepted tolerance for the generated coefficients.
///
/// Only the generator's main TU includes this (it is host-only: it spawns threads). The `-O0` net
/// builder units are compiled `-fno-exceptions -fno-rtti` and must not see it.
#pragma once

#include "numtracer/numeric/dpoly.hpp"
#include "numtracer/core/envvar.hpp" // env_flag / env_int — the single truth test for NT_* switches
#include "numtracer/numeric/mpoly.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <system_error>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace numtracer::numeric
{

  /// @brief Construct an empty polynomial of the phase-B backend type `P` at symbol-space size `nsym`.
  ///        For the real backends `MPoly`/`DPoly` this routes through the private-factory attorneys (the
  ///        `nsym`-taking ctors are not public); any OTHER `P` (e.g. a unit test's mock polynomial with a
  ///        public `explicit P(int)`) falls back to the plain ctor, so the fold templates stay generic.
  template <class P> inline P zero_like(int nsym)
  {
    if constexpr (std::is_same_v<P, MPoly>)
      return MPolyFactory::zero(nsym);
    else if constexpr (std::is_same_v<P, DPoly>)
      return DPolyFactory::zero(nsym);
    else
      return P(nsym);
  }

  /// @brief Scale a contracted trace by its sub-term scalar. Overloaded so phase B is one template
  ///        over both backends: the plain path multiplies by a constant `MPoly`, the dressed path
  ///        scales every kinematic coefficient (@ref scaleCx).
  inline MPoly scale_trace(int nsym, const MPoly &p, Cx c) { return MPolyFactory::scaled(nsym, p, c); }
  inline DPoly scale_trace(int, const DPoly &p, Cx c) { return scaleCx(p, c); }

  /// @brief Approximate heap footprint of a polynomial, for the trace-table RAM report. The table is
  ///        the one place this design trades memory for time, so the generator prints what it costs.
  inline std::size_t poly_bytes(const MPoly &p) { return p.terms.size() * sizeof(std::pair<Mono, Cx>); }
  inline std::size_t poly_bytes(const DPoly &p)
  {
    std::size_t b = 0;
    for (const auto &[d, m] : p.terms)
      b += d.size() * sizeof(int) + poly_bytes(m);
    return b;
  }

  /// @brief Run `f(i)` for `i` in `[0, n)` across `W` threads, dispatched off one flat atomic counter.
  ///        Work-stealing rather than a static split: even a uniform work list has a long tail, and
  ///        this is also what lets phase A ignore the per-net skew entirely. Returns the worker count
  ///        actually used. Falls back to fewer threads (down to the caller's own) if the system
  ///        refuses to spawn them.
  template <class F> unsigned parallel_flat(long n, unsigned W, F &&f)
  {
    if (n <= 0) return 0u;
    const unsigned nw = static_cast<unsigned>(std::min<long>(std::max(1u, W), n));
    std::atomic<long> next{0};
    auto work = [&] {
      long i;
      while ((i = next.fetch_add(1)) < n) f(i);
    };
    std::vector<std::thread> pool;
    pool.reserve(nw - 1);
    for (unsigned w = 1; w < nw; ++w) {
      try {
        pool.emplace_back(work);
      } catch (const std::system_error &) {
        break; // out of threads: the remaining work still runs, just on fewer of them
      }
    }
    work();
    for (auto &t : pool)
      t.join();
    return nw;
  }

  /// @brief PHASE A — contract traces `[0, nCache)` once each, in parallel over a flat work list.
  ///
  /// `nCache` is how many of the distinct traces are cached. Codegen.m orders the traces by
  /// *descending reference count*, so a cap keeps the traces that repay caching most; traces at index
  /// `>= nCache` are recomputed on demand in phase B (see @ref fold_net). Two reasons to cap:
  ///  - a trace referenced exactly once costs the same either way, so caching it is pure RAM for no
  ///    saving — Codegen.m's default `nCache` therefore excludes the singletons;
  ///  - `NT_GEN_MEMO_MAX` trims it further when memory is tight (the RAM lever: the dense flows are
  ///    memory-bound before they are compute-bound).
  ///
#if NT_PHASEA_STATS
  /// Print the merged phase-A counter/timer report plus the per-trace wall-time distribution.
  /// Stats builds only (see stats.hpp); called by contract_traces after the workers have joined, so
  /// the counters hold exactly phase A.
  inline void phasea_stats_report(const std::vector<double> &tsec, unsigned W)
  {
    const stats::PhaseACounters m = stats::merged();
    double tot = 0, mx = 0;
    std::size_t mxk = 0;
    for (std::size_t k = 0; k < tsec.size(); ++k) {
      tot += tsec[k];
      if (tsec[k] > mx) { mx = tsec[k]; mxk = k; }
    }
    // log10 histogram of per-trace wall times, buckets 10^-6..10^2 s
    int hist[9] = {0};
    for (double t : tsec) {
      int b = 0;
      for (double edge = 1e-6; b < 8 && t >= edge; ++b, edge *= 10.0) {}
      ++hist[b];
    }
    // top-10 traces by wall time
    std::vector<std::size_t> idx(tsec.size());
    for (std::size_t k = 0; k < idx.size(); ++k) idx[k] = k;
    const std::size_t K = std::min<std::size_t>(10, idx.size());
    std::partial_sort(idx.begin(), idx.begin() + K, idx.end(),
                      [&](std::size_t a, std::size_t b) { return tsec[a] > tsec[b]; });
    double topK = 0;
    for (std::size_t i = 0; i < K; ++i) topK += tsec[idx[i]];

    const double ts = 1e-9; // ns → s
    const double sect = ts * (m.t_dirac + m.t_elem + m.t_contract + m.t_reduce + m.t_divmono + m.t_divpoly);
    auto pct = [&](stats::nt_u64 ns_) { return sect > 0 ? 100.0 * ts * ns_ / sect : 0.0; };
    std::fprintf(stderr, "[stats] traces %llu  cpu %.1f s  max %.3f s (trace %zu)  top%zu %.1f s (%.1f%%)  makespan_lb(W=%u) %.1f s\n",
                 (unsigned long long)m.traces, tot, mx, mxk, K, topK, tot > 0 ? 100.0 * topK / tot : 0.0, W,
                 std::max(mx, W ? tot / W : tot));
    std::fprintf(stderr, "[stats] wall hist (log10 s): <1u %d  1-10u %d  10-100u %d  0.1-1m %d  1-10m %d  10-100m %d  0.1-1s %d  1-10s %d  >10s %d\n",
                 hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7], hist[8]);
    std::fprintf(stderr, "[stats] split: dirac %.1f s (%.1f%%)  elem %.1f s (%.1f%%)  contract %.1f s (%.1f%%)  reduce %.1f s (%.1f%%)  divmono %.1f s (%.1f%%)  divpoly %.1f s (%.1f%%)\n",
                 ts * m.t_dirac, pct(m.t_dirac), ts * m.t_elem, pct(m.t_elem), ts * m.t_contract, pct(m.t_contract),
                 ts * m.t_reduce, pct(m.t_reduce), ts * m.t_divmono, pct(m.t_divmono), ts * m.t_divpoly, pct(m.t_divpoly));
    std::fprintf(stderr, "[stats] contract detail: score %.1f s (%.1f%%)  in-step reduce %.1f s (%.1f%%)  elim arith %.1f s (%.1f%%)  steps %llu\n",
                 ts * m.t_cf_score, pct(m.t_cf_score), ts * m.t_cf_reduce, pct(m.t_cf_reduce),
                 ts * (m.t_contract - std::min(m.t_contract, m.t_cf_score + m.t_cf_reduce)),
                 pct(m.t_contract - std::min(m.t_contract, m.t_cf_score + m.t_cf_reduce)), (unsigned long long)m.cf_steps);
    std::fprintf(stderr, "[stats] mpoly: mul %llu (empty %llu, blocked %llu, prod-terms %llu)  add %llu (empty-copied %llu, empty-moved %llu)  from_scratch %llu (terms-in %llu)\n",
                 (unsigned long long)m.mul_calls, (unsigned long long)m.mul_empty, (unsigned long long)m.mul_blocked,
                 (unsigned long long)m.mul_prod_terms, (unsigned long long)m.add_calls, (unsigned long long)m.add_empty,
                 (unsigned long long)m.add_moved, (unsigned long long)m.fs_calls, (unsigned long long)m.fs_terms_in);
    std::fprintf(stderr, "[stats] reduce_units %llu calls (noop %llu) / %llu work-items  divmono %llu (noop %llu)  divpoly %llu (prefiltered %llu, trials %llu, exact %llu)\n",
                 (unsigned long long)m.ru_calls, (unsigned long long)m.ru_noop, (unsigned long long)m.ru_work,
                 (unsigned long long)m.dma_calls, (unsigned long long)m.dma_noop, (unsigned long long)m.dpa_calls,
                 (unsigned long long)m.dpa_pref, (unsigned long long)m.dpa_trials, (unsigned long long)m.dpa_exact);
    std::fprintf(stderr, "[stats] dirac: loops %llu (odd-skip %llu)  tokens %llu  4^f assignments %llu  mul2 %llu\n",
                 (unsigned long long)m.nd_calls, (unsigned long long)m.nd_odd_skip, (unsigned long long)m.nd_tokens,
                 (unsigned long long)m.nd_assign, (unsigned long long)m.mul2_calls);
  }
#endif // NT_PHASEA_STATS

  /// @param trace `trace(k) -> P`, the contraction of distinct trace `k`. Must be pure and safe to
  ///        call concurrently (the numeric_value_* entry points take everything by const reference).
  template <class P, class TraceFn>
  std::vector<P> contract_traces(int nsym, long nCache, unsigned W, TraceFn &&trace)
  {
    std::vector<P> T(static_cast<std::size_t>(std::max(0L, nCache)), zero_like<P>(nsym));
#if NT_PHASEA_STATS
    std::vector<double> tsec(T.size(), 0.0);
    parallel_flat(nCache, W, [&](long k) {
      const auto t0 = std::chrono::steady_clock::now();
      T[static_cast<std::size_t>(k)] = trace(static_cast<int>(k));
      tsec[static_cast<std::size_t>(k)] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    });
    phasea_stats_report(tsec, W);
#else
    parallel_flat(nCache, W, [&](long k) { T[static_cast<std::size_t>(k)] = trace(static_cast<int>(k)); });
#endif
    return T;
  }

  /// @brief PHASE B — fold one net: `Σ_j sc[j] · trace(idx[j])`.
  ///
  /// Reduced as a **binary-counter tree**: partial sums are kept on a stack tagged with a rank (a
  /// power of two = how many leaves it covers), and a new term carries into the stack exactly as a
  /// binary increment does, merging only equal ranks. So every addition combines two operands of
  /// comparable size, and at most O(log n) partials are ever live — where a left fold against one
  /// growing accumulator pays O(|acc|) per term, which on the heavy nets (>1000 terms) is the
  /// difference between O(n log n) and O(n²) element touches. Exactly `n - 1` additions either way.
  ///
  /// Traces at `idx[j] >= nCache` are not resident and are recomputed here; by construction those are
  /// referenced once, so nothing is computed twice.
  template <class P, class TraceFn>
  P fold_net(int nsym, const std::vector<int> &traceIdx, const std::vector<Cx> &subScale,
             const std::vector<P> &traceTable, long nCache, TraceFn &&trace)
  {
    std::vector<P> stack;
    std::vector<std::size_t> stackRank;
    P recomputed;

    for (std::size_t j = 0; j < traceIdx.size(); ++j) {
      const int k = traceIdx[j];
      // Bind, never copy: a ternary over `const P&` and a prvalue would materialise a copy of the
      // cached polynomial on every use — and the whole point is that these are used 5-8x each.
      const P *src;
      if (k < static_cast<int>(nCache)) {
        src = &traceTable[static_cast<std::size_t>(k)];
      } else {
        recomputed = trace(k);
        src = &recomputed;
      }

      P cur = scale_trace(nsym, *src, subScale[j]);
      std::size_t curRank = 1;
      while (!stack.empty() && stackRank.back() == curRank) {
        cur = stack.back() + cur; // earlier + later: the term order of the original left fold
        stack.pop_back();
        stackRank.pop_back();
        curRank *= 2;
      }
      stack.push_back(std::move(cur));
      stackRank.push_back(curRank);
    }

    if (stack.empty()) return zero_like<P>(nsym);
    P acc = std::move(stack.front()); // O(log n) leftovers, largest first
    for (std::size_t i = 1; i < stack.size(); ++i)
      acc = acc + stack[i];
    return acc;
  }

  /// @brief PHASE B (lever (b), dressed) — fold one net whose traces are PLAIN @ref MPoly into a
  ///        @ref DPoly: `Σ_j subScale[j] · traceTable[traceIdx[j]] ⊗ subDress[j]`.
  ///
  /// This is the dressed analogue of @ref fold_net. The trace table is plain `MPoly` (the dressing
  /// dimension was stripped at codegen time and lives in the per-sub-term `subDress`/`subScale`), so the
  /// same net can reference one concrete trace across MANY dressing channels without re-contracting it —
  /// the whole point of lever (b). Each sub-term becomes a one-term `DPoly` (`subDress[j]` →
  /// `subScale[j]·trace`), and those are summed by the SAME balanced binary-counter tree as @ref fold_net
  /// (so heavy nets stay O(n log n)); `DPoly::operator+` collects the channels. Traces at
  /// `traceIdx[j] >= nCache` are recomputed here, exactly as in @ref fold_net.
  template <class TraceFn>
  DPoly fold_net_dressed(int nsym, const std::vector<int> &traceIdx, const std::vector<Cx> &subScale,
                         const std::vector<DMono> &subDress, const std::vector<MPoly> &traceTable, long nCache,
                         TraceFn &&trace)
  {
    std::vector<DPoly> stack;
    std::vector<std::size_t> stackRank;
    MPoly recomputed;

    for (std::size_t j = 0; j < traceIdx.size(); ++j) {
      const int k = traceIdx[j];
      const MPoly *src;
      if (k < static_cast<int>(nCache)) {
        src = &traceTable[static_cast<std::size_t>(k)];
      } else {
        recomputed = trace(k);
        src = &recomputed;
      }

      // one dressing channel for this sub-term: subDress[j] · (subScale[j] · trace). Empty if the scaled
      // trace cancels to nothing (DPoly::add drops an empty MPoly), matching fold_net's zero handling.
      DPoly cur = DPolyFactory::zero(nsym);
      cur.add(subDress[j], scale_trace(nsym, *src, subScale[j]));
      std::size_t curRank = 1;
      while (!stack.empty() && stackRank.back() == curRank) {
        cur = stack.back() + cur; // earlier + later: preserve the original left-fold term order
        stack.pop_back();
        stackRank.pop_back();
        curRank *= 2;
      }
      stack.push_back(std::move(cur));
      stackRank.push_back(curRank);
    }

    if (stack.empty()) return zero_like<DPoly>(nsym);
    DPoly acc = std::move(stack.front());
    for (std::size_t i = 1; i < stack.size(); ++i)
      acc = acc + stack[i];
    return acc;
  }

  /// @brief PHASE B, driver — fold every net, in parallel over the nets.
  ///
  /// TEST-ONLY. Codegen.m no longer emits a call to this: production goes through
  /// @ref fold_groups_streaming, which folds each group on demand and drains straight to the sink so
  /// no net polynomial outlives its group. This one materialises them all, which is exactly why it
  /// makes a good reference — `tests/test_trace_fold.cpp` diffs the streaming fold against it and
  /// requires bit-identity. Keep it, but do not mistake it for a live path.
  ///
  /// Still net-parallel (unlike phase A), but that is fine here: the contraction is done, so a net's
  /// fold is proportional to its term count rather than to thousands of matrix products, and the
  /// merge in Codegen.m has already collapsed each net's repeated traces into one scalar apiece.
  template <class P, class TraceFn>
  std::vector<P> fold_nets(int nsym, const std::vector<std::vector<int>> &traceIdx,
                           const std::vector<std::vector<Cx>> &subScale, const std::vector<P> &traceTable,
                           long nCache, unsigned W, TraceFn &&trace)
  {
    std::vector<P> mp(traceIdx.size(), zero_like<P>(nsym));
    parallel_flat(static_cast<long>(traceIdx.size()), W, [&](long i) {
      const auto netIdx = static_cast<std::size_t>(i);
      mp[netIdx] = fold_net<P>(nsym, traceIdx[netIdx], subScale[netIdx], traceTable, nCache, trace);
    });
    return mp;
  }

  /// @brief How many NET polynomials may be in flight — and therefore held — at once.
  ///        `NT_GEN_GROUP_WINDOW` overrides; 0/unset picks `max(64, 8*W)`.
  ///
  /// This is the RAM lever for phase B, and it is a monotone dial: a window >= the net count lifts the
  /// bound entirely, reproducing the old all-resident schedule, so the pre-streaming behaviour stays
  /// reachable for A/B.
  ///
  /// The budget is counted in NETS, not groups, for two reasons. Nets are the unit that actually costs
  /// memory (a group accumulator is roughly one net's worth, while a group spans 1..30 nets depending
  /// on the flow), so nets are the honest currency for a RAM bound. And nets are the unit of
  /// PARALLELISM: @ref fold_groups_streaming dispatches a flat work list of nets, never of groups.
  /// Making the group the work item looks tempting — it needs no per-net storage — but it serialises
  /// every net inside a group onto one thread, which on a flow with few fat groups is catastrophic
  /// (measured on ZA4: 54 nets in 6 groups, one of them 27 nets, took the generator run 0.78 s -> 3.8 s).
  /// The flat net-level work list is the whole reason phase A scales — see this file's header note on
  /// per-net skew — and phase B must not give it up.
  inline long net_window(long nNet, unsigned W)
  {
    static const long ov = env_int("NT_GEN_GROUP_WINDOW", 0); // 0 = no override, use the rule below
    if (nNet <= 0) return 0;
    if (ov > 0) return std::min(ov, nNet);
    return std::min(nNet, std::max<long>(64, 8L * static_cast<long>(W)));
  }

  /// @brief Does `groups` cover every net exactly once? Pure predicate, so it is unit-testable
  ///        without running a generator; @ref check_group_partition is the reporting wrapper.
  ///
  /// @ref fold_groups_streaming folds each group's members on demand, so it is only equivalent to
  /// `fold_nets` + an eager group loop if `groups` PARTITIONS the nets — which Codegen.m guarantees by
  /// construction (`Complement`/`GatherBy` over all net indices). A duplicate would mean silently
  /// folding a net twice; a gap would mean silently dropping one from the kernel. Both are cheap to
  /// detect here and expensive to debug downstream.
  inline bool group_partition_stats(const std::vector<std::vector<int>> &groups, long nNet, long &dup,
                                    long &gap, long &oob)
  {
    std::vector<int> cov(static_cast<std::size_t>(std::max(0L, nNet)), 0);
    dup = gap = oob = 0;
    for (const auto &grp : groups)
      for (int d : grp) {
        if (d < 0 || d >= nNet) ++oob;
        else ++cov[static_cast<std::size_t>(d)];
      }
    for (int c : cov) {
      if (c > 1) ++dup;
      else if (c == 0) ++gap;
    }
    return dup == 0 && gap == 0 && oob == 0;
  }

  /// @brief Abort the generation if `groups` does not partition the nets.
  ///
  /// FATAL, and run on EVERY generation. It used to be emitted behind `if(ntprof)` — i.e. gated on
  /// NT_GEN_PROFILE, which no production regeneration sets — and it only warned to stderr, where the
  /// line was lost in the generator log. So the one guard standing between a mis-built group list and
  /// a silently wrong kernel had, in practice, never run at all. It is O(nNet) once per generation
  /// (nets are tens to low thousands) against a build measured in seconds to minutes.
  ///
  /// Verified before being made fatal: a full regeneration of all 29 flows in DEFAULT_FLOWS reported
  /// zero violations, so no committed flow legitimately produces a non-partition.
  ///
  /// Exits rather than throwing: a nonzero exit is already the generator's wired-up failure signal
  /// (`MakeNTKernel::genfail`, which tests/gen/regen_check.sh greps for), and parts of the generator
  /// are compiled `-fno-exceptions`.
  inline void check_group_partition(const std::vector<std::vector<int>> &groups, long nNet)
  {
    long dup = 0, gap = 0, oob = 0;
    if (group_partition_stats(groups, nNet, dup, gap, oob)) return;
    std::fprintf(stderr,
                 "[num] FATAL groups is not a partition: %ld duplicated, %ld uncovered, %ld "
                 "out-of-range (of %ld nets) — a duplicate folds a net into the kernel twice, a gap "
                 "drops one; either is a silently wrong kernel\n",
                 dup, gap, oob, nNet);
    std::exit(1);
  }

  /// @brief PHASE B, streaming driver — fold nets in bounded waves and hand each group's accumulator
  ///        straight to `sink`, so no net polynomial outlives the group that consumes it.
  ///
  /// `fold_nets` materialises one fully-expanded polynomial per net and returns them all; the caller
  /// then sums them into per-group accumulators. Both sets are live at once and neither is ever
  /// released, which on the dense 4-point flows is 20+ GB (488 nets x ~41 MB) against a 20 MB trace
  /// table — a ~1000x ratio between what is *cached* and what is merely *retained*. Nothing is
  /// revisited: each net polynomial is written once, read once by its group, then dead. So this is a
  /// streaming workload, and running it as a batch is the entire bug.
  ///
  /// Here at most `window` net polynomials are ever live, regardless of the flow's shape.
  ///
  /// Equivalence with `fold_nets` + an eager left fold is exact, not approximate:
  ///  - each group still sums its members left-to-right in group order over the same `fold_net`
  ///    results, so the floating-point association is unchanged bit-for-bit;
  ///  - `sink` is called on the CALLING thread for `gi = 0, 1, 2, ...` strictly ascending, exactly once
  ///    per group, which is what keeps `GlobalEnv` symbol-intern order and the shared CSE builder's
  ///    instruction stream identical to the non-streamed lowering.
  ///
  /// @param scale `scale(d, P&&) -> P` applies net `d`'s colour weight. Deliberately a caller-supplied
  ///        callable rather than @ref scale_trace: the emitted MPoly branches write
  ///        `mp[d] * env.constant(colv[d])` (poly x constant) while `scale_trace` computes
  ///        constant x poly. `MPoly::operator*` sort-collects so the two almost certainly agree — but a
  ///        refactor that must not move a single ulp should not rest on "almost certainly", and a
  ///        lambda preserves each branch's exact expression while keeping `colv` out of this signature.
  /// @param trace as in @ref fold_net: pure and safe to call concurrently.
  ///
  // Forward declaration: the two drivers below are thin wrappers over this shared skeleton (defined
  // just after them).
  template <class P, class FoldScaledFn, class Sink>
  void fold_groups_streaming_impl(int nsym, const std::vector<std::vector<int>> &groups, unsigned W,
                                  long window, FoldScaledFn &&foldScaledNet, Sink &&sink);

  /// The wave/streaming SKELETON is factored into @ref fold_groups_streaming_impl so the plain and the
  /// dressed (lever (b)) drivers share it: the only thing that varies between them is how one net's
  /// polynomial is produced (`fold_net<P>` + colour scale vs `fold_net_dressed` + colour scale), which
  /// each passes as a `foldScaledNet(netId) -> P` callable. The wrapper below reproduces the previous
  /// body's operations exactly, so the non-dressed emitted kernel is byte-identical.
  template <class P, class TraceFn, class ScaleFn, class Sink>
  void fold_groups_streaming(int nsym, const std::vector<std::vector<int>> &traceIdx,
                             const std::vector<std::vector<Cx>> &subScale,
                             const std::vector<std::vector<int>> &groups, const std::vector<P> &traceTable,
                             long nCache, unsigned W, long window, TraceFn &&trace, ScaleFn &&scale,
                             Sink &&sink)
  {
    fold_groups_streaming_impl<P>(
        nsym, groups, W, window,
        [&](int netId) {
          const auto netIdx = static_cast<std::size_t>(netId);
          return scale(netId,
                       fold_net<P>(nsym, traceIdx[netIdx], subScale[netIdx], traceTable, nCache, trace));
        },
        std::forward<Sink>(sink));
  }

  /// @brief PHASE B, streaming driver — lever (b) dressed variant. Trace table is plain @ref MPoly; each
  ///        net folds into a @ref DPoly via @ref fold_net_dressed (carrying the per-sub-term dressing
  ///        monomials @p subDress), then the colour @p scale and @p sink run exactly as in the plain
  ///        driver. Same bounded live set and ascending-group sink order, so `GlobalEnv` intern order /
  ///        the CSE instruction stream are governed identically.
  template <class TraceFn, class ScaleFn, class Sink>
  void fold_groups_streaming_dressed(int nsym, const std::vector<std::vector<int>> &traceIdx,
                                     const std::vector<std::vector<Cx>> &subScale,
                                     const std::vector<std::vector<DMono>> &subDress,
                                     const std::vector<std::vector<int>> &groups,
                                     const std::vector<MPoly> &traceTable, long nCache, unsigned W,
                                     long window, TraceFn &&trace, ScaleFn &&scale, Sink &&sink)
  {
    fold_groups_streaming_impl<DPoly>(
        nsym, groups, W, window,
        [&](int netId) {
          const auto netIdx = static_cast<std::size_t>(netId);
          return scale(netId, fold_net_dressed(nsym, traceIdx[netIdx], subScale[netIdx], subDress[netIdx],
                                               traceTable, nCache, trace));
        },
        std::forward<Sink>(sink));
  }

  /// @brief The shared wave/streaming skeleton of the phase-B drivers. @p foldScaledNet returns the
  ///        colour-scaled polynomial of the net it is given; @p sink absorbs each group's accumulator in
  ///        ascending group order on the calling thread. Bounded live set: at most `window` net
  ///        polynomials in flight.
  template <class P, class FoldScaledFn, class Sink>
  void fold_groups_streaming_impl(int nsym, const std::vector<std::vector<int>> &groups, unsigned W,
                                  long window, FoldScaledFn &&foldScaledNet, Sink &&sink)
  {
    const long nG = static_cast<long>(groups.size());
    if (nG <= 0) return;

    // Per-wave RSS trace (NT_GEN_PROFILE=2). Phase B's peak has three possible shapes and they call for
    // opposite fixes, so measure rather than guess: a MONOTONE climb means something accumulates across
    // waves (the lowered programs, the shared CSE env, or allocator churn) and the window is irrelevant;
    // a SAWTOOTH with a flat ceiling means the peak is the transient of the heaviest trace contraction
    // in a wave, which only a lower worker count can reduce; a climb that tracks wave size means the
    // window itself is the cost.
    const bool wprof = env_int("NT_GEN_PROFILE", 0) == 2;
    // Any-level profile: the phase-B wall is really three unlike costs — the parallel per-net fold,
    // the serial group-sum drain, and the serial lowering hiding inside `sink`. One fused number
    // cannot say which of them a slow flow is paying, so split them here.
    const bool pprof = env_flag("NT_GEN_PROFILE");
    double foldSec = 0.0, drainSec = 0.0, sinkSec = 0.0;
    using pclock = std::chrono::steady_clock;
    auto secsSince = [](pclock::time_point t0) { return std::chrono::duration<double>(pclock::now() - t0).count(); };
    auto rssMB = [] {
      long tot = 0, res = 0;
      if (FILE *f = std::fopen("/proc/self/statm", "r")) {
        if (std::fscanf(f, "%ld %ld", &tot, &res) != 2) res = 0;
        std::fclose(f);
      }
      return res * 4096.0 / 1048576.0; // page size is 4 KiB where this runs
    };
    long waveNo = 0;

    // Work in WAVES of consecutive groups, each sized by a budget of NETS rather than of groups, so a
    // flow with a few fat groups and one with many thin ones both get a bounded live set and the same
    // flat work list. A wave always takes at least one group, so an oversized group still makes
    // progress — its nets are simply all in flight at once, which is the irreducible floor.
    std::vector<P> part;   // per-net folded+scaled results for this wave, in wave-local order
    std::vector<long> off; // off[k] = index into `part` where wave-local group k starts
    std::vector<int> flat; // this wave's nets, flattened — the parallel work list

    for (long g0 = 0; g0 < nG;) {
      long g1 = g0, nets = 0;
      while (g1 < nG) {
        const long sz = static_cast<long>(groups[static_cast<std::size_t>(g1)].size());
        if (g1 > g0 && window > 0 && nets + sz > window) break;
        nets += sz;
        ++g1;
      }

      flat.clear();
      off.clear();
      off.push_back(0);
      for (long gi = g0; gi < g1; ++gi) {
        for (int d : groups[static_cast<std::size_t>(gi)])
          flat.push_back(d);
        off.push_back(static_cast<long>(flat.size()));
      }

      // FLAT NET-LEVEL PARALLELISM — the same granularity `fold_nets` used, and the reason this
      // scales: sub-terms per net are wildly skewed, so scheduling per group (or anything coarser than
      // a net) lets one fat work item pin utilisation.
      part.clear();
      part.reserve(flat.size());
      for (std::size_t i = 0; i < flat.size(); ++i)
        part.push_back(zero_like<P>(nsym));
      const double rssPre = wprof ? rssMB() : 0.0;
      const auto tFold = pclock::now();
      parallel_flat(static_cast<long>(flat.size()), W, [&](long i) {
        const auto j = static_cast<std::size_t>(i);
        part[j] = foldScaledNet(flat[j]);
      });
      foldSec += secsSince(tFold);
      const double rssFold = wprof ? rssMB() : 0.0;

      // Drain the wave in ascending group order, on THIS thread.
      for (long gi = g0; gi < g1; ++gi) {
        const long k = gi - g0;
        const auto tDrain = pclock::now();
        P acc = zero_like<P>(nsym);
        for (long i = off[static_cast<std::size_t>(k)]; i < off[static_cast<std::size_t>(k + 1)]; ++i) {
          const auto j = static_cast<std::size_t>(i);
          acc = acc + part[j];
          part[j] = zero_like<P>(nsym); // release as absorbed
        }
        drainSec += secsSince(tDrain);
        const auto tSink = pclock::now();
        sink(static_cast<std::size_t>(gi), std::move(acc));
        sinkSec += secsSince(tSink);
      }

      if (wprof) {
        // pre -> after-fold -> after-sink. The pre->fold delta is what the wave's nets plus the
        // concurrent trace contractions cost; the fold->sink delta is what lowering retained; and the
        // drift of `pre` across waves is what never comes back.
        std::fprintf(stderr, "[num]   wave %ld: groups %ld..%ld, %zu nets, RSS %.0f -> %.0f -> %.0f MB\n",
                     waveNo, g0, g1 - 1, flat.size(), rssPre, rssFold, rssMB());
        std::fflush(stderr);
      }
      ++waveNo;

      g0 = g1;
    }

    if (pprof) {
      std::fprintf(stderr,
                   "[num]   phase B split: net-fold %.1f s (W=%u), group-sum %.1f s (serial), lower %.1f s (serial)\n",
                   foldSec, W, drainSec, sinkSec);
      std::fflush(stderr);
    }
  }

} // namespace numtracer::numeric
