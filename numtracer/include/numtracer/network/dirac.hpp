/// @file dirac.hpp
/// @brief The closed gamma-chain tokens (@ref DFac / @ref DiracChain) the generator builds, plus
///        `dirac_value`, a Wick-pairing trace used as a TEST ORACLE.
///
/// The production trace is `numeric_dirac` (numeric/numeric_contract.hpp): it multiplies the chain
/// as 4×4 Weyl-block matrices of polynomials and handles every token kind. `dirac_value` instead
/// folds a FREE/SLASH-only chain into a Lorentz @ref LorentzNet (metrics / vectors over the free legs)
/// by the Wick pairing recursion — an independent algorithm that tests/test_numeric_contract.cpp
/// grades `numeric_dirac` against. It refuses every token kind it cannot trace.
///
/// A chain is a list of trace-ordered tokens (already cyclically closed by the front-end's
/// `orderDiracFacs`, mathematica/CodegenNets.m): each is either a FREE gluon leg `γ^μ` (an open
/// Lorentz id `mu`, contracts the projector later), a SLASHED propagator `γ·p` (the momentum
/// `p = Σ coeff·fund(vid)` as a `vlc`, mirroring @ref LorentzFactor's vector linear combination), or a γ5
/// marker.
///
/// The trace is the standard Wick pairing recursion
///   tr(t1 … t_{2n}) = Σ_{j≥2} (−1)^j  g(t1,t_j)  tr(t1 … t̂_j … t_{2n}),
/// with `g(free μa, free μb)=δ` → @ref metric, `g(free μ, slash p)=p^μ` → @ref vec, and
/// `g(slash p, slash q)=p·q` emitted as two @ref vec leaves on a **fresh** internal Lorentz label so
/// the numeric contraction sums the two-vector index class into one scalar product. The overall `tr(1)=4` is
/// folded into the coefficient. An odd number of gammas traces to 0.
#pragma once

#include "numtracer/network/network.hpp"

#include <utility>
#include <vector>

namespace numtracer::inline network
{

  /// @brief One token of a closed, trace-ordered gamma chain (a @ref DiracChain).
  ///   - `gamma(mu)`  : a free γ^μ — `mu` is an OPEN Lorentz index, contracted by the Lorentz network
  ///   - `slash(k)`   : a slashed momentum γ·k
  ///   - `gamma5()`   : γ5 (chiral trace → Levi-Civita)
  ///   - `chargeC()`  : the charge-conjugation matrix C
  ///   - `comm(A, B)` : the BARE commutator `[A,B] = A·B − B·A` of two gammas as ONE chain token.
  ///                    Each leg is a free γ (open index) or a slashed momentum. It keeps the
  ///                    antisymmetric γ-pair folded so the commutator is never distributed into two
  ///                    separate traces. The σ^{μν}=(i/2)[γ^μ,γ^ν] normalization (the i/2 and any
  ///                    sign) lives in the scalar coefficient, not the token.
  ///   - `loop_sep()` : separates two independent spinor loops in one chain
  ///   - `transposed(t)` : the transpose of token `t`
  // Members are non-const so DFac stays movable in the std::vector chains it lives in (const members
  // would delete move-assignment and force the vlc vectors to be copied on every reallocation). The
  // builder functions below are the only constructors, so the values are still effectively immutable.
  struct DFac {
    enum Kind { Gamma, Gamma5, Slash, Comm, LoopSep, C };
    Kind kind = Gamma;
    // Field docs name the VARIANT, never its enum ordinal: ordinals shift whenever `Kind` grows.
    int mu = -1; ///< Gamma: open Lorentz id. Comm: leg-A FREE id (-1 ⇒ leg-A is a slash, use `vlc`)
    Vlc vlc;     ///< Slash: the momentum lin. comb. Comm: leg-A slash momentum (when mu < 0)
    int nu = -1; ///< Comm: leg-B FREE id (-1 ⇒ leg-B is a slash, use `vlc2`)
    /// @brief Multiply this factor's TRANSPOSE instead of the factor itself.
    ///
    /// A closed spinor loop is a cycle in a degree-2 index network, and following that cycle is a
    /// matrix product only where each factor's declared (din,dout) order matches the traversal.
    /// Where it does not, the network wants @f$M^T@f$ — that is not a gamma identity, it is what
    /// following the cycle means. The front-end walk (`orderDiracFacs`) knows which factors it
    /// entered backwards and marks them; this flag carries that marking. Handled entirely at
    /// BLOCK-PRECOMPUTE time in `numeric_dirac`, so the fold itself is unaffected.
    ///
    /// Placed HERE, between `nu` and `vlc2`, so it lands in existing padding: sizeof(DFac) stays 64
    /// (72 at the end cost +1.4% instructions in the fold).
    bool transposed = false;
    Vlc vlc2; ///< Comm: leg-B slash momentum (when nu < 0)
  };
  /// @brief A closed gamma chain in trace order (loop closes implicitly).
  using DiracChain = std::vector<DFac>;

  /// @brief A free γ^μ: its Lorentz index stays open until the Lorentz network contracts it.
  inline DFac gamma(LorentzIndex mu) { return {DFac::Gamma, mu.id, {}, -1, false, {}}; }
  /// @brief Boundary between two INDEPENDENT closed spinor loops in one chain (e.g. a quark loop
  ///        and the projection-closed external line, tied together only by gluon propagators). Each
  ///        loop is traced separately and the resulting Lorentz tensors are contracted through the
  ///        Lorentz network — NOT a Wick pairing across loops.
  inline DFac loop_sep() { return {DFac::LoopSep, -1, {}, -1, false, {}}; }
  /// @brief A slashed momentum `k̸ = γ·k`.
  inline DFac slash(Momentum k) { return {DFac::Slash, -1, std::move(k.lc), -1, false, {}}; }
  /// @brief γ5 (chiral trace → Levi-Civita).
  inline DFac gamma5() { return {DFac::Gamma5, -1, {}, -1, false, {}}; }
  /// @brief The charge-conjugation matrix @f$C=\gamma^2\gamma^4@f$ as a chain token.
  ///
  /// Like γ5 it is block-DIAGONAL in the Weyl split, so it does NOT flip the antidiagonal trace
  /// parity — but unlike γ5 its blocks are signed permutations rather than ±identity, so it costs
  /// two 2×2 multiplies rather than a sign flip. It carries no index of its own: the two spinor
  /// slots are the chain neighbours, exactly as for γ5.
  ///
  /// A `C` reaching the engine means the front end chose to KEEP it rather than fold it away; see
  /// the charge-conjugation rewrite in `CodegenNets.m`. Both are legal — the folded form is production,
  /// the token form is the oracle they are graded against.
  inline DFac chargeC() { return {DFac::C, -1, {}, -1, false, {}}; }

  /// @brief The transpose of a chain token: `transposed(gamma(mu))`, `transposed(slash(k))`, …
  ///
  /// Composes with every builder above rather than doubling them. In the Weyl block split a
  /// transpose is cheap and exact: a block-ANTIdiagonal factor (γ, slash) transposes by swapping its
  /// P and Q blocks and transposing each 2×2; a block-DIAGONAL one (γ5, C, σ) transposes its two
  /// diagonal blocks in place. Special cases the engine exploits: @f$\gamma_5^T=\gamma_5@f$ is a
  /// no-op and @f$C^T=-C@f$ is a sign flip.
  inline DFac transposed(DFac d)
  {
    d.transposed = !d.transposed;
    return d;
  }
  /// @brief The bare commutator `[A, B] = AB − BA` of two gammas, kept as ONE token. Each leg is a
  ///        free `γ^μ` (a @ref LorentzIndex) or a slashed momentum (a @ref Momentum). The
  ///        σ^{μν} = (i/2)[γ^μ,γ^ν] normalisation is NOT included.
  inline DFac comm(LorentzIndex mu, LorentzIndex nu) { return {DFac::Comm, mu.id, {}, nu.id, false, {}}; }
  inline DFac comm(Momentum a, Momentum b) { return {DFac::Comm, -1, std::move(a.lc), -1, false, std::move(b.lc)}; }
  inline DFac comm(LorentzIndex mu, Momentum b) { return {DFac::Comm, mu.id, {}, -1, false, std::move(b.lc)}; }
  inline DFac comm(Momentum a, LorentzIndex nu) { return {DFac::Comm, -1, std::move(a.lc), nu.id, false, {}}; }

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
      if (af && bf) return metric(LorentzIndex{a.mu}, LorentzIndex{b.mu});    // δ^{μa μb}
      if (af && !bf) return vec_lc(a.mu, b.vlc); // p_b^{μa}
      if (!af && bf) return vec_lc(b.mu, a.vlc); // p_a^{μb}
      const int sharedLbl = fresh++;             // slash–slash → p_a · p_b (sum over a fresh shared index)
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
      if (tokens.empty()) return {LorentzTerm{Cx{1, 0}, {}}}; // empty trace: tr(1)/4 = 1 (tr(1)=4 applied by dirac_value)
      LorentzNet acc;                                       // empty == 0
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
    // rather than return a plausible-looking wrong answer. (The production engine, numeric_dirac,
    // handles both.)
    for (const DFac &d : chain)
      if (d.kind == DFac::Comm)
        NT_THROW(std::runtime_error,
                 "dirac_value: chain contains a Comm (sigma) token, which the Wick-pairing "
                 "recursion cannot trace (pair_factor would read its empty vlc and silently "
                 "collapse the trace). Use numeric_dirac for sigma-bearing chains.");
      else if (d.kind == DFac::LoopSep)
        NT_THROW(std::runtime_error,
                 "dirac_value: chain contains a LoopSep marker — split the chain into its "
                 "independent spinor loops before tracing (see ndetail::split_loops).");
      else if (d.transposed)
        NT_THROW(std::runtime_error,
                 "dirac_value: chain contains a TRANSPOSED token. The Wick-pairing recursion pairs "
                 "factors by the Clifford algebra and has no notion of a factor's index order, so it "
                 "would silently trace the untransposed factor. Use numeric_dirac.");
      else if (d.kind == DFac::Gamma5)
        NT_THROW(std::runtime_error,
                 "dirac_value: chain contains a gamma5 token. The Wick-pairing recursion knows only "
                 "the four Clifford generators — pair_factor would read gamma5's empty vlc as a "
                 "zero slash and silently return a wrong trace. Use numeric_dirac.");
      else if (d.kind == DFac::C)
        NT_THROW(std::runtime_error,
                 "dirac_value: chain contains a charge-conjugation (C) token. The Wick-pairing "
                 "recursion knows only the Clifford algebra, and C is not a gamma — pair_factor "
                 "would read its empty vlc and silently collapse the trace. Use numeric_dirac, "
                 "or fold C away in the front end before tracing here.");
    // Gamma parity. Count ONLY the antidiagonal-block tokens, exactly as numeric_dirac does: a Comm
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

} // namespace numtracer::network
