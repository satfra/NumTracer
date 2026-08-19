/// @file env.hpp
/// @brief @ref numtracer::numeric::LorentzEnv — a factory that binds the symbol-space size `nsym`
///        (and the unit-constraint groups) ONCE and is the sole sanctioned construction path for
///        @ref MPoly / @ref DPoly.
///
/// Every polynomial fed into one trace must share an `nsym` (the length of each monomial's exponent
/// vector); mismatches were previously guarded only by debug asserts in the arithmetic operators and
/// were easy to introduce (a bare integer threaded through dozens of call sites). `LorentzEnv` captures
/// `nsym` once so a wrong value is unconstructible at the API surface: `env.var(i)`, `env.constant(c)`,
/// `env.numeric_value_netval(...)` all bake in the env's `nsym`. It is a construction-time factory only
/// — the polynomials it mints still carry a bare inline `int nsym`, so the arithmetic hot path is
/// unchanged (no env pointer is stored per object). It is a `friend` of @ref MPoly / @ref DPoly, so it
/// reaches their (private) `nsym`-taking factories directly.
///
/// The colour/flavour analogue is @ref numtracer::network::SUNEnv (binds the group rank), in
/// `network/sun_net.hpp`.
#pragma once

#include "numtracer/numeric/numeric_contract.hpp" // MPoly/DPoly/Mat4, NNet, numeric_value* entry points
#include "numtracer/numeric/numeric_driver.hpp"   // mpoly_to_cpp

#include <array>
#include <utility>
#include <vector>

namespace numtracer::numeric
{

  // Forward declarations of the host-only trace-fold templates (`numeric/trace_fold.hpp`). That header
  // spawns threads and is forbidden in the -fno-exceptions net-builder TUs, so `env.hpp` must NOT include
  // it. The two @ref LorentzEnv methods that forward here are member templates whose bodies are only
  // instantiated in the generator's main TU, where `trace_fold.hpp` is already included — so these
  // declarations are enough to compile `env.hpp` everywhere else.
  template <class P, class TraceFn>
  std::vector<P> contract_traces(int nsym, long nCache, unsigned W, TraceFn &&trace);
  template <class P, class TraceFn>
  std::vector<P> fold_nets(int nsym, const std::vector<std::vector<int>> &traceIdx,
                           const std::vector<std::vector<Cx>> &subScale, const std::vector<P> &traceTable,
                           long nCache, unsigned W, TraceFn &&trace);
  template <class P, class TraceFn, class ScaleFn, class Sink>
  void fold_groups_streaming(int nsym, const std::vector<std::vector<int>> &traceIdx,
                             const std::vector<std::vector<Cx>> &subScale,
                             const std::vector<std::vector<int>> &groups, const std::vector<P> &traceTable,
                             long nCache, unsigned W, long window, TraceFn &&trace, ScaleFn &&scale,
                             Sink &&sink);
  template <class TraceFn, class ScaleFn, class Sink>
  void fold_groups_streaming_dressed(int nsym, const std::vector<std::vector<int>> &traceIdx,
                                     const std::vector<std::vector<Cx>> &subScale,
                                     const std::vector<std::vector<DMono>> &subDress,
                                     const std::vector<std::vector<int>> &groups,
                                     const std::vector<MPoly> &traceTable, long nCache, unsigned W,
                                     long window, TraceFn &&trace, ScaleFn &&scale, Sink &&sink);

  /// @brief Binds `nsym` (+ unit groups) once; the sole sanctioned construction path for @ref MPoly /
  ///        @ref DPoly and the env-aware form of the numeric backend's public entry points.
  ///
  /// Cheap value object (an int + a small vector), so multiple symbol spaces coexist naturally: a
  /// diagram that needs a second frame just constructs a second `LorentzEnv`.
  class LorentzEnv
  {
  public:
    explicit LorentzEnv(int nsym, std::vector<std::vector<int>> units = {})
        : nsym_(nsym), units_(std::move(units))
    {
    }

    int nsym() const { return nsym_; }

    // ── polynomial construction (forwards to the private MPoly/DPoly factories; friend access) ────────
    MPoly zero() const { return MPoly(nsym_); }
    MPoly constant(Cx c) const { return MPoly::constant(nsym_, c); }
    MPoly var(int i) const { return MPoly::var(nsym_, i); }
    MPoly mono(const std::vector<int> &e, Cx c) const { return MPoly::mono(nsym_, e, c); }
    MPoly atom(int aid) const { return MPoly::atom(nsym_, aid); }
    DPoly dzero() const { return DPoly(nsym_); }

    // ── spinor-matrix layer (`spinor_mat.hpp`) ───────────────────────────────────────────────────────
    // Fully qualified because this member name shadows the free function of the same name.
    Mat4 slashC(const std::array<MPoly, 4> &comp) const { return ::numtracer::numeric::slashC(nsym_, comp); }

    // ── public entry points (env supplies nsym_, and units_ to the *_netval forms) ───────────────────
    MPoly numeric_value(const network::DiracNet &dirac, const NNet &lorentz,
                        const std::vector<std::array<MPoly, 4>> &comp, const std::vector<MPoly> &atomDen) const
    {
      return ::numtracer::numeric::numeric_value(nsym_, dirac, lorentz, comp, atomDen);
    }
    MPoly numeric_value_netval(const network::DiracNet &dirac, const network::NetVal &lor,
                               const std::vector<std::array<MPoly, 4>> &comp,
                               const std::vector<MPoly> &atomDen) const
    {
      return ::numtracer::numeric::numeric_value_netval(nsym_, dirac, lor, comp, atomDen, units_);
    }
    DPoly numeric_value_dressed(const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                const NNet &lorentz, const std::vector<std::array<MPoly, 4>> &comp,
                                const std::vector<MPoly> &atomDen) const
    {
      return ::numtracer::numeric::numeric_value_dressed(nsym_, chain, slots, lorentz, comp, atomDen);
    }
    DPoly numeric_value_dressed_netval(const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                       const network::NetVal &lor, const std::vector<std::array<MPoly, 4>> &comp,
                                       const std::vector<MPoly> &atomDen) const
    {
      return ::numtracer::numeric::numeric_value_dressed_netval(nsym_, chain, slots, lor, comp, atomDen, units_);
    }
    MPoly numeric_value_dressed_netval_mp(const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                          const network::NetVal &lor, const std::vector<std::array<MPoly, 4>> &comp,
                                          const std::vector<MPoly> &atomDen) const
    {
      return ::numtracer::numeric::numeric_value_dressed_netval_mp(nsym_, chain, slots, lor, comp, atomDen, units_);
    }
    std::vector<MPoly> collect_atom_denoms(const std::vector<network::NetVal> &lors,
                                           const std::vector<std::array<MPoly, 4>> &comp) const
    {
      return ::numtracer::numeric::collect_atom_denoms(nsym_, lors, comp);
    }

    // ── trace-fold phases (member templates; free templates forward-declared above) ──────────────────
    template <class P, class TraceFn>
    std::vector<P> contract_traces(long nCache, unsigned W, TraceFn &&trace) const
    {
      return ::numtracer::numeric::contract_traces<P>(nsym_, nCache, W, std::forward<TraceFn>(trace));
    }
    template <class P, class TraceFn>
    std::vector<P> fold_nets(const std::vector<std::vector<int>> &traceIdx,
                             const std::vector<std::vector<Cx>> &subScale, const std::vector<P> &traceTable,
                             long nCache, unsigned W, TraceFn &&trace) const
    {
      return ::numtracer::numeric::fold_nets<P>(nsym_, traceIdx, subScale, traceTable, nCache, W,
                                                std::forward<TraceFn>(trace));
    }
    /// Streaming phase B: folds each group's nets on demand and drains straight to `sink`, so no net
    /// polynomial outlives its group. See `trace_fold.hpp` for the equivalence argument.
    template <class P, class TraceFn, class ScaleFn, class Sink>
    void fold_groups_streaming(const std::vector<std::vector<int>> &traceIdx,
                               const std::vector<std::vector<Cx>> &subScale,
                               const std::vector<std::vector<int>> &groups, const std::vector<P> &traceTable,
                               long nCache, unsigned W, long window, TraceFn &&trace, ScaleFn &&scale,
                               Sink &&sink) const
    {
      ::numtracer::numeric::fold_groups_streaming<P>(nsym_, traceIdx, subScale, groups, traceTable, nCache, W,
                                                     window, std::forward<TraceFn>(trace),
                                                     std::forward<ScaleFn>(scale), std::forward<Sink>(sink));
    }
    /// Streaming phase B, lever (b) dressed variant: a plain-MPoly `traceTable` + the per-sub-term
    /// dressing monomials `subDress` fold into a DPoly per net. See `trace_fold.hpp`'s
    /// `fold_groups_streaming_dressed`.
    template <class TraceFn, class ScaleFn, class Sink>
    void fold_groups_streaming_dressed(const std::vector<std::vector<int>> &traceIdx,
                                       const std::vector<std::vector<Cx>> &subScale,
                                       const std::vector<std::vector<DMono>> &subDress,
                                       const std::vector<std::vector<int>> &groups,
                                       const std::vector<MPoly> &traceTable, long nCache, unsigned W,
                                       long window, TraceFn &&trace, ScaleFn &&scale, Sink &&sink) const
    {
      ::numtracer::numeric::fold_groups_streaming_dressed(nsym_, traceIdx, subScale, subDress, groups,
                                                          traceTable, nCache, W, window,
                                                          std::forward<TraceFn>(trace),
                                                          std::forward<ScaleFn>(scale), std::forward<Sink>(sink));
    }

  private:
    int nsym_;                             ///< number of user symbols (MPoly variable count)
    std::vector<std::vector<int>> units_;  ///< unit-constraint groups (ΣUμ²=1), threaded to *_netval
  };

} // namespace numtracer::numeric
