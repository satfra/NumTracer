/// @file contract_dressed.hpp
/// @brief Contraction of diagrams with dressed numerators: sums of Dirac structures whose
///        coefficients are run-time dressing functions.
///
/// Dressed numerators: a vertex or propagator numerator may be a sum of Dirac structures, each times
/// a dressing function known only at kernel run time (e.g. A(p)·p̸ + B(p), or Σ_i λ_i(p)·T_i^μ). It
/// enters the chain as a SLOT (@ref DSlot, one @ref DSlotOpt per term). The slot combinations are
/// enumerated and each concrete chain goes through the same pipeline; @ref ndetail::contract_dressed
/// collects the results per dressing monomial (a @ref DPoly), while @ref ndetail::contract_structural
/// sums them into a plain @ref Poly for the generator, which carries the dressings separately
/// (see engine/generator_phases.hpp).
///
///
/// A dressed propagator numerator like `Mq·δ + Z(p)·(γ·p)` is a SUM of Dirac structures whose
/// per-structure coefficients are *runtime* dressings. The front-end keeps such a sum eager (it no
/// longer distributes the diagram into one copy per structure), and the generator hands it here as a
/// Dirac chain with SLOT tokens: each slot is the structure sum of one dressed numerator. We collect
/// it WITHOUT the `2^D`-diagram blowup: enumerate the structure choices, contract each concrete chain
/// with the EXISTING @ref ndetail::contract (so all the validated Dirac/Lorentz machinery is reused
/// verbatim), and accumulate the results into ONE @ref DPoly keyed by the dressing monomial. The whole
/// diagram then lowers to ONE trace function whose dressing factors the shared CSE/Horner collects.
///
/// Layout: declarations first (always visible), then the bodies, which compile only in the library
/// TU or a header-only build (see core/export.hpp).
#pragma once

#include "numtracer/core/export.hpp"     // NUMTRACER_FUNC / NUMTRACER_DEFINE_BODIES
#include "numtracer/engine/contract.hpp" // ndetail::contract_net, InputScan (bodies)
#include "numtracer/poly/dpoly.hpp"      // DPoly / DMono

#include <utility>
#include <vector>

namespace numtracer::inline engine
{

  // ================================================================================================
  // Declarations: slot types and the entry points (always visible)
  // ================================================================================================

  /// @brief One structure option of a collected Dirac slot: a numeric
  ///        coefficient × a product of dressing atoms (`dress`) × a Dirac structure, where the
  ///        structure is
  ///          - a **Dirac-token chain** @ref toks spliced in place of the slot (spinor `din→dout`);
  ///            its FREE-Lorentz tokens (a `gamma(μ)` / an open `comm` leg) are open legs on the
  ///            Dirac side, and
  ///          - a set of **Lorentz-net factors** @ref netFacs (a vector `p^μ`, a metric `g^{μν}`, …)
  ///            carrying any remaining open legs, appended to the surrounding net.
  ///
  /// The union of the two sets' free Lorentz ids is the slot's shared open-leg set `{μ₁,…,μ_k}`, closed
  /// by the surrounding net — the SAME set for every option (so the net contracts a fixed leg set
  /// regardless of structure choice). This one form covers every case, any leg count `k ≥ 0`:
  ///  - propagator numerator `δ` ⇒ `toks={}, netFacs={}`; slash `γ·p̸` ⇒ `toks={slash(vlc)}` (k=0);
  ///  - single-gluon vertex T1 `γ^μ` ⇒ `toks={gamma(μ)}`; T4 `p̸₁γ^μ` ⇒ `toks={slash(p1),gamma(μ)}`;
  ///    T7 `σ^{μν}p̸_ν` ⇒ `toks={comm(μ,vlc)}` (k=1);
  ///  - open leg on a vector `p^μ·δ` ⇒ `netFacs={Vector(μ,p)}`; two-gluon `g^{μν}·δ` ⇒
  ///    `netFacs={Metric(μ,ν)}`; two open γ's ⇒ `toks={gamma(μ),gamma(ν)}` (k=2), etc.
  struct DSlotOpt {
    Cx coeff{1, 0};
    std::vector<int> dress;             ///< dressing-atom ids (need not be pre-sorted; merged sorted on use)
    std::vector<network::DFac> toks;    ///< Dirac-token chain spliced for this option (free-Lorentz tokens = open legs)
    std::vector<LorentzFactor> netFacs; ///< extra Lorentz-net factors for this option (open legs on vectors/metrics/…)
  };
  /// @brief A dressed numerator = the sum of its structure options.
  using DSlot = std::vector<DSlotOpt>;
  /// @brief One token of a dressed Dirac chain: either a FIXED @ref network::DFac (γ / slash / γ5 /
  ///        commutator) or a SLOT reference (an index into the diagram's slot list).
  struct DChainTok {
    bool isSlot = false;
    network::DFac fac{network::DFac::Gamma, -1, {}, -1, false, {}};
    int slot = -1;
    /// @brief Splice this slot's chain TRANSPOSED (slot tokens only; a fixed token carries its own
    ///        `DFac::transposed` instead).
    ///
    /// A slot is a sum of chains, and @f$(A_1\cdots A_n)^T = A_n^T\cdots A_1^T@f$ — so transposing
    /// one means reversing each option's token list and transposing every token in it. That has to
    /// happen HERE, at the splice, rather than in the front end: `$dsCache`/`$dslCache` key on the
    /// slot expression, so a Mathematica-side transposed slot would collide with its untransposed
    /// twin in the cache. The option's `netFacs` are pure Lorentz tensors on the open legs and are
    /// unaffected by a spinor transpose.
    bool transposed = false;
  };
  inline DChainTok dtfix(network::DFac f) { return {false, std::move(f), -1, false}; }
  inline DChainTok dtslot(int s) { return {true, {network::DFac::Gamma, -1, {}, -1, false, {}}, s, false}; }
  /// @brief A slot reference spliced transposed. See @ref DChainTok::transposed.
  inline DChainTok dtrslot(int s) { return {true, {network::DFac::Gamma, -1, {}, -1, false, {}}, s, true}; }

  namespace ndetail
  {
    /// @brief Dressed-slot analogue of @ref contract: returns the collected @ref DPoly (one channel
    ///        per dressing monomial). The engine behind the dressed @ref Frame::contract.
    NUMTRACER_FUNC DPoly contract_dressed(int nsym, const std::vector<DChainTok> &chain,
                                          const std::vector<DSlot> &slots, const LorentzNet &lor,
                                          const std::vector<std::array<Poly, 4>> &comp,
                                          const std::vector<Poly> &atomDen,
                                          const std::vector<std::vector<int>> &units = {});
    /// @brief STRUCTURAL variant of @ref contract_dressed: contract the collected slots into a PLAIN
    ///        @ref Poly, DISCARDING the dressing dimension (the dressing monomial keys are summed
    ///        away). The generator strips each option's dressing (`coeff`→1, `dress`→{}) into a
    ///        per-sub-term scalar/monomial and feeds the dressing-free structure here, so the trace
    ///        table dedups on structure alone (dressing variants of one concrete trace collapse) and
    ///        each entry is a plain `Poly` with no dressing dimension. The 4^(#collapsed loop)
    ///        tr(1)=4 factor is STRUCTURAL (it depends only on the concrete chain) and is kept in the
    ///        returned trace. Only meaningful when every option's `dress` is empty (the generator
    ///        guarantees this); if fed genuine dressings it would silently sum them.
    NUMTRACER_FUNC Poly contract_structural(int nsym, const std::vector<DChainTok> &chain,
                                            const std::vector<DSlot> &slots, const LorentzNet &lor,
                                            const std::vector<std::array<Poly, 4>> &comp,
                                            const std::vector<Poly> &atomDen,
                                            const std::vector<std::vector<int>> &units = {});
  } // namespace ndetail

#if NUMTRACER_DEFINE_BODIES

  // ================================================================================================
  // 1. Enumerating slot combinations
  // ================================================================================================

  /// @brief Append a collected-slot combination's open-leg net factors (@ref DSlotOpt::netFacs — a
  ///        vector `p^μ`, a metric `g^{μν}`, …) to EVERY term of a Lorentz net, so the surrounding net
  ///        closes their legs. This is the exact Lorentz structure the distributed diagram emits, which
  ///        is why it can reuse all the validated contraction machinery unchanged.
  ///
  /// An empty @p facs returns @p lor untouched — the case for every option whose structure is purely
  /// Dirac-side (the fast path).
  ///
  /// An empty @p lor is SEEDED with a unit product term first. "Append to every term" is a no-op on a
  /// net with no terms, so appending straight into `lor` would DROP @p facs — silently, since the
  /// result is still empty and the contraction then takes its `lor.empty()` shortcut into
  /// @ref close_free_legs, where the option's now-unclosed legs abort (or, when they happen to pair
  /// up, self-contract to a wrong number). An empty net means "the Lorentz rest is the scalar 1", and
  /// that is exactly a default-constructed `LorentzTerm` (coefficient 1, no factors) — the same
  /// thing the front end emits as `konst(1.0)`.
  inline LorentzNet with_slot_facs(const LorentzNet &lor, const std::vector<LorentzFactor> &facs)
  {
    if (facs.empty()) return lor;
    LorentzNet out = lor.empty() ? LorentzNet{LorentzTerm{}} : lor;
    for (LorentzTerm &pt : out)
      pt.e.insert(pt.e.end(), facs.begin(), facs.end());
    return out;
  }

  namespace ndetail
  {
    /// Enumerate the Cartesian product of slot-structure choices, contract each concrete Dirac chain
    /// via @p contract (an `Poly`-returning closure that runs the full Dirac+Lorentz contraction), and
    /// collect the results into one @ref DPoly keyed by the dressing monomial. Each slot is referenced
    /// at most once in @p chain (a numerator occupies one chain position). @p contract receives the
    /// concrete Dirac chain AND this combination's extra Lorentz-net factors (the chosen options'
    /// @ref DSlotOpt::netFacs — open legs on vectors/metrics/…); empty for every option whose structure
    /// is purely Dirac-side ⇒ the caller takes a byte-identical fast path.
    ///
    /// The enumeration is factored out of @ref dress_collect so the two accumulation policies share it:
    /// the DPoly collection (`emit(dmono, mp)` → `out.add(dmono, mp)`) and the STRUCTURAL Poly reduction
    /// (@ref dress_collect_mp: `emit(_, mp)` → `out += mp`, dropping the dressing key). @p emit is called
    /// once per surviving combination with the sorted dressing monomial and the (coeff·trace) `Poly`.
    template <class ContractFn, class Emit>
    inline void dress_enumerate(int nsym, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                ContractFn &&contract, Emit &&emit)
    {
      const int nSlots = static_cast<int>(slots.size());
      // The dressed chain is a quark line with ≥1 closed spinor loop ((LoopSep markers)+1). A combination
      // that collapses a WHOLE loop to the identity (every slot on it takes the δ option, no fixed γ) leaves
      // an empty segment: split_loops drops it and close_loops returns 1, silently losing that loop's
      // tr(1)=4 — a 4× undercount (e.g. the all-mass channel of a Yukawa self-energy quark loop). Count the
      // loops up front so each combination can restore 4^(#collapsed loops). Flows whose every loop keeps a
      // fixed γ (Zq/ZAqbq{1,4,7}, …) never collapse (nEmpty==0).
      int nloops = 1;
      for (const DChainTok &tok : chain)
        if (!tok.isSlot && tok.fac.kind == network::DFac::LoopSep) ++nloops;
      std::vector<int> choice(nSlots, 0);
      bool done = false;
      while (!done) {
        // build the concrete chain + accumulate this combination's numeric coeff and dressing monomial
        Cx combCoeff{1, 0};
        DMono dressMono;
        network::DiracChain concrete;
        concrete.reserve(chain.size());
        std::vector<LorentzFactor> extraNet; // this combination's options' net factors (usually empty)
        for (const DChainTok &tok : chain) {
          if (!tok.isSlot) {
            concrete.push_back(tok.fac);
            continue;
          }
          const DSlotOpt &opt = slots[tok.slot][choice[tok.slot]];
          combCoeff = combCoeff * opt.coeff;
          dressMono.insert(dressMono.end(), opt.dress.begin(), opt.dress.end());
          // General collected Dirac slot: splice this option's Dirac-token chain in place of
          // the slot (its free-Lorentz tokens are open legs), and collect its Lorentz-net factors — the
          // surrounding net closes every open leg for all structure choices alike. An option with an
          // empty `toks` is the spinor identity δ (tr(1)=4 restored by nCollapsed below when a whole
          // loop collapses); empty `netFacs` ⇒ no net change (the fast path).
          if (tok.transposed) {
            // (A1..An)^T = An^T..A1^T: reverse the option's chain and transpose every token in it.
            // netFacs are pure Lorentz and untouched. Loop accounting below is order-insensitive
            // (it counts LoopSep markers and split_loops segments), so reversing is safe here.
            for (auto it = opt.toks.rbegin(); it != opt.toks.rend(); ++it)
              concrete.push_back(network::transposed(*it));
          } else {
            concrete.insert(concrete.end(), opt.toks.begin(), opt.toks.end());
          }
          extraNet.insert(extraNet.end(), opt.netFacs.begin(), opt.netFacs.end());
        }
        if (!(combCoeff.re == 0 && combCoeff.im == 0)) {
          // restore tr(1)=4 for every spinor loop that collapsed to the identity in this combination
          const int nCollapsed = nloops - static_cast<int>(split_loops(concrete).size());
          for (int c = 0; c < nCollapsed; ++c)
            combCoeff = combCoeff * Cx{4, 0};
          Poly mp = contract(concrete, extraNet);
          if (!mp.empty()) {
            mp = mp * PolyFactory::constant(nsym, combCoeff);
            emit(dmono_sorted(std::move(dressMono)), std::move(mp));
          }
        }
        // odometer over the slot choices: advance the lowest slot, carrying into the next on wrap.
        if (nSlots == 0) break;
        int k = 0;
        for (; k < nSlots; ++k) {
          if (++choice[k] < static_cast<int>(slots[k].size())) break;
          choice[k] = 0;
        }
        if (k == nSlots) done = true;
      }
    }

    /// @brief The DPoly collection: accumulate each combination into its dressing-monomial channel.
    template <class ContractFn>
    inline DPoly dress_collect(int nsym, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                               ContractFn &&contract)
    {
      DPoly out = DPolyFactory::zero(nsym);
      dress_enumerate(nsym, chain, slots, std::forward<ContractFn>(contract),
                      [&](const DMono &d, Poly &&mp) { out.add(d, mp); });
      return out;
    }

    /// @brief The STRUCTURAL Poly reduction: sum every combination into ONE plain Poly,
    ///        discarding the dressing monomial. Correct only when the slots carry no dressing (the
    ///        generator strips it out first); see @ref contract_structural. The debug
    ///        assert below makes that precondition checkable instead of comment-only: fed genuine
    ///        dressings, this would silently sum structures that belong in different channels.
    ///
    /// The left-to-right accumulation is deliberate and is NOT a quadratic hazard: the generator
    /// expands each structure×dressing combination into its own SINGLE-option sub-term
    /// (mathematica/CodegenGenerator.m), so every caller has exactly one combination. A tree fold would
    /// reassociate the sums (≤ 1 ulp) and could shift the emitted literals for no gain.
    template <class ContractFn>
    inline Poly dress_collect_mp(int nsym, const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                 ContractFn &&contract)
    {
#ifndef NDEBUG
      for (const DSlot &s : slots)
        for (const DSlotOpt &o : s)
          assert(o.dress.empty() && "dress_collect_mp requires dressing-free slots (lever (b))");
#endif
      Poly out = PolyFactory::zero(nsym);
      dress_enumerate(nsym, chain, slots, std::forward<ContractFn>(contract),
                      [&](const DMono &, Poly &&mp) { out = std::move(out) + std::move(mp); });
      return out;
    }
  } // namespace ndetail

  // ================================================================================================
  // 2. Entry points: contract_dressed, contract_structural
  // ================================================================================================

  namespace ndetail
  {
    /// @ref InputScan over a dressed diagram: the fixed tokens, EVERY slot option (so each structure
    /// combination is covered by one scan) and the Lorentz net.
    inline void check_dressed_inputs(const std::vector<DChainTok> &chain, const std::vector<DSlot> &slots,
                                     const LorentzNet &lor, std::size_t nMom, std::size_t nAtom)
    {
      InputScan scan;
      for (const DChainTok &tok : chain)
        if (!tok.isSlot) scan.add(tok.fac);
      for (const DSlot &slot : slots)
        for (const DSlotOpt &opt : slot) {
          for (const DFac &d : opt.toks)
            scan.add(d);
          for (const LorentzFactor &f : opt.netFacs)
            scan.add(f);
        }
      scan.add(lor);
      scan.verify(nMom, nAtom);
    }

    NUMTRACER_FUNC DPoly contract_dressed(int nsym, const std::vector<DChainTok> &chain,
                                          const std::vector<DSlot> &slots, const LorentzNet &lor,
                                          const std::vector<std::array<Poly, 4>> &comp,
                                          const std::vector<Poly> &atomDen, const std::vector<std::vector<int>> &units)
    {
      check_dressed_inputs(chain, slots, lor, comp.size(), atomDen.size());
      return dress_collect(nsym, chain, slots, [&](const DiracChain &d, const std::vector<LorentzFactor> &slotFacs) {
        return slotFacs.empty() ? contract_net(nsym, d, lor, comp, atomDen, units)
                                : contract_net(nsym, d, with_slot_facs(lor, slotFacs), comp, atomDen, units);
      });
    }

    NUMTRACER_FUNC Poly contract_structural(int nsym, const std::vector<DChainTok> &chain,
                                            const std::vector<DSlot> &slots, const LorentzNet &lor,
                                            const std::vector<std::array<Poly, 4>> &comp,
                                            const std::vector<Poly> &atomDen,
                                            const std::vector<std::vector<int>> &units)
    {
      check_dressed_inputs(chain, slots, lor, comp.size(), atomDen.size());
      return dress_collect_mp(nsym, chain, slots, [&](const DiracChain &d, const std::vector<LorentzFactor> &slotFacs) {
        return slotFacs.empty() ? contract_net(nsym, d, lor, comp, atomDen, units)
                                : contract_net(nsym, d, with_slot_facs(lor, slotFacs), comp, atomDen, units);
      });
    }
  } // namespace ndetail

#endif // NUMTRACER_DEFINE_BODIES

} // namespace numtracer::inline engine
