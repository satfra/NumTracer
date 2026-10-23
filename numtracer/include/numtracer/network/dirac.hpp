/// @file dirac.hpp
/// @brief The closed gamma-chain tokens (@ref DFac / @ref DiracChain) the generator builds.
///
/// A chain is a list of trace-ordered tokens (already cyclically closed by the front-end's
/// `orderDiracFacs`, mathematica/CodegenNets.m): a FREE gluon leg `γ^μ` (an open Lorentz id `mu`,
/// contracted by the Lorentz network later), a SLASHED propagator `γ·p` (the momentum
/// `p = Σ coeff·fund(vid)` as a `vlc`, mirroring @ref LorentzFactor's vector linear combination), γ5,
/// the charge-conjugation matrix C, a bare commutator `[A,B]`, or a separator between two independent
/// spinor loops. The trace itself is `trace_dirac_loop` (engine/contract.hpp), which multiplies the
/// chain as 4×4 Weyl-block matrices of polynomials.
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
    /// BLOCK-PRECOMPUTE time in `trace_dirac_loop`, so the fold itself is unaffected.
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

} // namespace numtracer::inline network
