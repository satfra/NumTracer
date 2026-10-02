/// @file stats.hpp
/// @brief Phase-A instrumentation: opt-in counters + section timers for the contraction engine.
///
/// Compile-time gated on `NT_PHASEA_STATS`. When the macro is absent (the default, and the only
/// configuration ever installed/shipped) every hook expands to `((void)0)` and this header adds no
/// includes beyond nothing — the engine object code is byte-identical. A stats build is a
/// throwaway profiling binary: compile EVERY TU of that binary with `-DNT_PHASEA_STATS=1`
/// (mixing stats and non-stats TUs is an ODR violation on the inline hot paths — the same failure
/// mode as the documented stale-libNumTracer.a bug) and never install it.
///
/// Include constraints (same as the rest of the numeric engine): this header is pulled into
/// `mpoly.hpp`/`numeric_contract.hpp`, which the `-O0 -fno-exceptions` net-builder TUs may see —
/// so no `<thread>`, no `<mutex>`, no exceptions. Aggregation uses an atomic intrusive list of
/// heap-allocated (deliberately leaked) per-thread counter blocks; the dump walks the list after
/// phase A's workers have joined (trace_fold.hpp owns the dump).
#pragma once

#if NT_PHASEA_STATS

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>

namespace numtracer::numeric::stats
{

  using nt_u64 = std::uint64_t;

  /// One thread's phase-A counters. Fields are totals; times are nanoseconds.
  struct PhaseACounters {
    // MPoly arithmetic (mpoly.hpp)
    nt_u64 mul_calls = 0;      ///< operator* invocations with both sides non-empty
    nt_u64 mul_empty = 0;      ///< operator* short-circuits (an operand was the zero poly)
    nt_u64 mul_prod_terms = 0; ///< Σ |a|·|b| product monomials emitted into scratch
    nt_u64 mul_blocked = 0;    ///< operator* calls that took the chunked (>kMulMaxScratch) path
    nt_u64 mul_keyed = 0;      ///< operator* calls on atom-free inline operands (key-only scratch)
    nt_u64 mul_const = 0;      ///< operator* calls with a single-constant factor (no scratch, no sort)
    nt_u64 add_calls = 0;      ///< operator+ invocations with both sides non-empty
    nt_u64 add_empty = 0;      ///< operator+ short-circuits that COPY the surviving side (lvalue)
    nt_u64 add_moved = 0;      ///< operator+ short-circuits that MOVED the surviving side (rvalue)
    nt_u64 ru_noop = 0;        ///< reduce_units pass-through early-exits (nothing to rewrite)
    nt_u64 dma_noop = 0;       ///< divThroughMonomialAtoms pass-through early-exits (no atoms)
    nt_u64 dpa_pref = 0;       ///< divThroughPolyAtoms trials rejected by the lead pre-filter (no map copy)
    nt_u64 fs_calls = 0;       ///< from_scratch invocations
    nt_u64 fs_terms_in = 0;    ///< Σ scratch terms entering from_scratch (sort size)
    nt_u64 ru_calls = 0;       ///< reduce_units invocations (non-empty groups)
    nt_u64 ru_work = 0;        ///< Σ work-stack items popped in reduce_units (branch blowup)
    nt_u64 dma_calls = 0;      ///< divThroughMonomialAtoms invocations
    nt_u64 dpa_calls = 0;      ///< divThroughPolyAtoms invocations
    nt_u64 dpa_trials = 0;     ///< trial divisions attempted (each pays a Grp map copy today)
    nt_u64 dpa_exact = 0;      ///< trial divisions that were exact
    nt_u64 dpa_noop = 0;       ///< divThroughPolyAtoms pass-throughs (no multi-term atom present)
    // Dirac fold (numeric_contract.hpp)
    nt_u64 nd_calls = 0;       ///< numeric_dirac invocations (spinor loops)
    nt_u64 nd_tokens = 0;      ///< Σ chain tokens
    nt_u64 nd_assign = 0;      ///< Σ 4^f free-leg assignments folded
    nt_u64 nd_odd_skip = 0;    ///< odd-parity chains skipped (zero tensor, no arithmetic)
    nt_u64 mul2_calls = 0;     ///< 2×2 Weyl-block multiplies (8 MPoly mults + 4 adds each)
    nt_u64 traces = 0;         ///< numeric_value_netval invocations
    // Section times inside numeric_value_netval (ns)
    nt_u64 t_dirac = 0;    ///< dirac_loop_factors (the 4^f chain folds)
    nt_u64 t_elem = 0;     ///< elem_to_nelem + fuse_projectors + elem_factor building
    nt_u64 t_contract = 0; ///< contract_factors (greedy elimination incl. its in-step reductions)
    nt_u64 t_cf_score = 0; ///< within t_contract: the min-fill id scoring rescans
    nt_u64 t_cf_reduce = 0; ///< within t_contract: eliminate's per-outFlat reduce_units+divmono
    nt_u64 cf_steps = 0;   ///< elimination steps taken (one id summed out per step)
    nt_u64 t_reduce = 0;   ///< final reduce_units on the result
    nt_u64 t_divmono = 0;  ///< final divThroughMonomialAtoms
    nt_u64 t_divpoly = 0;  ///< final divThroughPolyAtoms
    PhaseACounters *next = nullptr; ///< intrusive registry link

    void add(const PhaseACounters &o)
    {
      const nt_u64 *s = &o.mul_calls;
      nt_u64 *d = &mul_calls;
      // PhaseACounters is standard-layout with uint64 fields up to `next`; sum them field-wise.
      constexpr std::size_t n = offsetof(PhaseACounters, next) / sizeof(nt_u64);
      for (std::size_t i = 0; i < n; ++i)
        d[i] += s[i];
    }
  };

  inline std::atomic<PhaseACounters *> &registry()
  {
    static std::atomic<PhaseACounters *> head{nullptr};
    return head;
  }

  /// This thread's counter block: heap-allocated on first use and pushed onto the registry, and
  /// deliberately leaked so the block outlives the worker thread (the dump runs after join).
  inline PhaseACounters &tls()
  {
    thread_local PhaseACounters *c = nullptr;
    if (!c) {
      c = new PhaseACounters();
      PhaseACounters *h = registry().load(std::memory_order_relaxed);
      do {
        c->next = h;
      } while (!registry().compare_exchange_weak(h, c, std::memory_order_release, std::memory_order_relaxed));
    }
    return *c;
  }

  /// Merge every registered thread's counters (call after phase-A workers joined).
  inline PhaseACounters merged()
  {
    PhaseACounters m;
    for (PhaseACounters *c = registry().load(std::memory_order_acquire); c; c = c->next)
      m.add(*c);
    return m;
  }

  /// RAII section timer accumulating into one `tls()` field.
  class Timer
  {
    nt_u64 PhaseACounters::*field_;
    std::chrono::steady_clock::time_point t0_;

  public:
    explicit Timer(nt_u64 PhaseACounters::*field) : field_(field), t0_(std::chrono::steady_clock::now()) {}
    ~Timer()
    {
      tls().*field_ += static_cast<nt_u64>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t0_).count());
    }
  };

} // namespace numtracer::numeric::stats

#define NT_STAT_ADD(field, v) (::numtracer::numeric::stats::tls().field += (v))
#define NT_STAT_TIMER(field) ::numtracer::numeric::stats::Timer nt_stat_timer_##field(&::numtracer::numeric::stats::PhaseACounters::field)

#else // !NT_PHASEA_STATS

#define NT_STAT_ADD(field, v) ((void)0)
#define NT_STAT_TIMER(field) ((void)0)

#endif // NT_PHASEA_STATS
