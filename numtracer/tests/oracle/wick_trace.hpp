/// @file wick_trace.hpp
/// @brief Test-only Wick-pairing Dirac trace (`dirac_value`), an independent oracle for the
///        production trace `trace_dirac_loop` (engine/contract.hpp).
///
/// `dirac_value` folds a FREE/SLASH-only chain into a Lorentz @ref numtracer::LorentzNet (metrics /
/// vectors over the free legs) by the standard Wick pairing recursion
///   tr(t1 … t_{2n}) = Σ_{j≥2} (−1)^j  g(t1,t_j)  tr(t1 … t̂_j … t_{2n}),
/// with `g(free μa, free μb)=δ` → `metric`, `g(free μ, slash p)=p^μ` → `vec`, and
/// `g(slash p, slash q)=p·q` emitted as two `vec` leaves on a **fresh** internal Lorentz label so the
/// numeric contraction sums the two-vector index class into one scalar product. The overall
/// `tr(1)=4` is folded into the coefficient; an odd number of gammas traces to 0. It refuses every
/// token kind it cannot trace (σ, γ5, C, loop separators, transposed tokens).
///
/// Graded against in tests/test_contract.cpp; its refusals are checked in tests/test_guards.cpp.
/// It is not part of the shipped library.
#pragma once

#include "numtracer/network/dirac.hpp"
#include "numtracer/network/network.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace numtracer::inline network
{

  namespace dirac_detail
  {

    /// @brief A vector leg carrying a full momentum linear combination `vlc` on Lorentz index `lbl`.
    inline LorentzNet vec_lc(int lbl, const std::vector<std::pair<double, int>> &vlc)
    {
      return {LorentzTerm{Cx{1, 0}, {LorentzFactor{.kind = LorentzFactor::Vector, .a = lbl, .b = -1, .vlc = vlc}}}};
    }

    /// @brief The pairing factor `g(a,b)` between two gamma tokens. Slash–slash uses a fresh shared
    ///        Lorentz label so the numeric contraction extracts the scalar product.
    inline LorentzNet pair_factor(const DFac &a, const DFac &b, int &fresh)
    {
      const bool af = (a.kind == DFac::Gamma), bf = (b.kind == DFac::Gamma);
      if (af && bf) return metric(LorentzIndex{a.mu}, LorentzIndex{b.mu}); // δ^{μa μb}
      if (af && !bf) return vec_lc(a.mu, b.vlc);                           // p_b^{μa}
      if (!af && bf) return vec_lc(b.mu, a.vlc);                           // p_a^{μb}
      const int sharedLbl = fresh++; // slash–slash → p_a · p_b (sum over a fresh shared index)
      return mul(vec_lc(sharedLbl, a.vlc), vec_lc(sharedLbl, b.vlc));
    }

    /// @brief Wick pairing recursion for a token list of FREE/SLASH gammas (no γ5). Returns the trace
    ///        divided by 4 (the `tr(1)=4` is applied by @ref dirac_value).
    ///
    /// Pin the first token and pair it with each later token j: tr = Σ_{j≥1} (−1)^{j+1} g(t0,tj)·tr(rest),
    /// where `rest` is the chain with t0 and tj removed and g(·,·) is @ref pair_factor. Base case: the
    /// empty chain traces to 1 (= tr(1)/4).
    inline LorentzNet trace_rec(const std::vector<DFac> &tokens, int &fresh)
    {
      if (tokens.empty())
        return {LorentzTerm{Cx{1, 0}, {}}}; // empty trace: tr(1)/4 = 1 (tr(1)=4 applied by dirac_value)
      LorentzNet acc;                       // empty == 0
      for (std::size_t j = 1; j < tokens.size(); ++j) {
        std::vector<DFac> rest;
        rest.reserve(tokens.size() - 2);
        for (std::size_t k = 1; k < tokens.size(); ++k)
          if (k != j) rest.push_back(tokens[k]);
        const double sgn = (j % 2 == 1) ? 1.0 : -1.0; // (−1)^{j+1}, 0-indexed (matches gammaTraceSum)
        LorentzNet pairFac = pair_factor(tokens[0], tokens[j], fresh);
        LorentzNet subTrace = trace_rec(rest, fresh);
        acc = add(std::move(acc), scale(sgn, mul(pairFac, subTrace)));
      }
      return acc;
    }

  } // namespace dirac_detail

  /// @brief Contract a closed gamma chain into a Lorentz @ref LorentzNet over its free legs.
  /// @param chain          the trace-ordered tokens (free legs / slashes only; anything else throws).
  /// @param firstFreeLabel a Lorentz id strictly above every label the surrounding component uses, so
  ///                       the fresh slash–slash labels never collide with the free legs / projector.
  /// @return the trace as a `LorentzNet` (empty == structural zero, e.g. an odd gamma count).
  inline LorentzNet dirac_value(const DiracChain &chain, int firstFreeLabel)
  {
    // trace_rec implements the Wick pairing for FREE/SLASH tokens only: pair_factor treats any
    // non-Gamma token as a slash and reads its `vlc`, which for a Comm (σ) built from two FREE legs
    // is EMPTY — vec_lc then folds to the zero 4-vector and the trace silently collapses. LoopSep is
    // not a gamma at all and would be paired as one. Neither can be traced here, so refuse them
    // rather than return a plausible-looking wrong answer. (The production engine, trace_dirac_loop,
    // handles both.)
    for (const DFac &d : chain)
      if (d.kind == DFac::Comm)
        NT_THROW(std::runtime_error, "dirac_value: chain contains a Comm (sigma) token, which the Wick-pairing "
                                     "recursion cannot trace (pair_factor would read its empty vlc and silently "
                                     "collapse the trace). Use trace_dirac_loop for sigma-bearing chains.");
      else if (d.kind == DFac::LoopSep)
        NT_THROW(std::runtime_error, "dirac_value: chain contains a LoopSep marker — split the chain into its "
                                     "independent spinor loops before tracing (see ndetail::split_loops).");
      else if (d.transposed)
        NT_THROW(std::runtime_error,
                 "dirac_value: chain contains a TRANSPOSED token. The Wick-pairing recursion pairs "
                 "factors by the Clifford algebra and has no notion of a factor's index order, so it "
                 "would silently trace the untransposed factor. Use trace_dirac_loop.");
      else if (d.kind == DFac::Gamma5)
        NT_THROW(std::runtime_error,
                 "dirac_value: chain contains a gamma5 token. The Wick-pairing recursion knows only "
                 "the four Clifford generators — pair_factor would read gamma5's empty vlc as a "
                 "zero slash and silently return a wrong trace. Use trace_dirac_loop.");
      else if (d.kind == DFac::C)
        NT_THROW(std::runtime_error, "dirac_value: chain contains a charge-conjugation (C) token. The Wick-pairing "
                                     "recursion knows only the Clifford algebra, and C is not a gamma — pair_factor "
                                     "would read its empty vlc and silently collapse the trace. Use trace_dirac_loop, "
                                     "or fold C away in the front end before tracing here.");
    // Gamma parity. Count ONLY the antidiagonal-block tokens, exactly as trace_dirac_loop does: a Comm
    // is two gammas, a LoopSep none and a C two, so all contribute 0 mod 2. Keep the two counts in
    // step even though those kinds are refused above.
    std::size_t nAntidiag = 0;
    for (const DFac &d : chain)
      if (d.kind == DFac::Gamma || d.kind == DFac::Slash) ++nAntidiag;
    if (nAntidiag % 2 == 1) return {}; // odd chain → 0
    int fresh = firstFreeLabel;
    LorentzNet r = dirac_detail::trace_rec(chain, fresh); // no γ5 tokens (refused above)
    return scale(Cx{4.0, 0}, std::move(r));               // tr(1) = 4
  }

} // namespace numtracer::inline network
