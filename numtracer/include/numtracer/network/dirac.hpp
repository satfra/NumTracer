/// @file dirac.hpp
/// @brief Numeric Dirac gamma-trace for the **build-time generator** — contracts a closed spinor
///        (gamma) loop into a Lorentz @ref NetVal (metrics / vectors / Levi-Civita over the free
///        gluon legs), at codegen time, so the Mathematica front-end no longer expands the trace.
///
/// This is the Dirac analogue of @ref sun_net.hpp's `SUNNet`: colour folds to a scalar `Cx`, but a
/// Dirac trace `tr(γ^{μ1} … γ^{μn})` is a **Lorentz tensor** over the free legs `μ_i` whose
/// coefficients are momentum-dependent (the slashed propagator momenta). So `dirac_value` folds to a
/// `NetVal` (not a scalar): the same metric/vector structure the symbolic `gammaTraceSum` produced
/// in Mathematica, but built in C++ and fed straight into the numeric contraction
/// (@ref numtracer::numeric::numeric_value_netval). Moving it here removes the
/// `(2n−1)!!`-terms-×-(hundreds-of-distributed-branches) explosion from the Mathematica codegen
/// (a 3-gluon-vertex quark triangle was minutes/GBs there).
///
/// A chain is a list of trace-ordered tokens (already cyclically closed by the front-end's
/// `orderDiracChain`): each is either a FREE gluon leg `γ^μ` (an open Lorentz id `mu`, contracts the
/// projector later), a SLASHED propagator `γ·p` (the momentum `p = Σ coeff·fund(vid)` as a `vlc`,
/// mirroring @ref Elem's vector linear combination), or a γ5 marker.
///
/// The trace is the standard Wick pairing recursion
///   tr(t1 … t_{2n}) = Σ_{j≥2} (−1)^j  g(t1,t_j)  tr(t1 … t̂_j … t_{2n}),
/// with `g(free μa, free μb)=δ` → @ref met, `g(free μ, slash p)=p^μ` → @ref vec, and
/// `g(slash p, slash q)=p·q` emitted as two @ref vec leaves on a **fresh** internal Lorentz label so
/// the numeric contraction sums the two-vector index class into one scalar product. The overall `tr(1)=4` is
/// folded into the coefficient. An odd number of gammas traces to 0. (Generator-only — runtime.)
#pragma once

#include "numtracer/network/network.hpp"

#include <utility>
#include <vector>

namespace numtracer::network
{

  /// @brief One token of a closed, trace-ordered gamma chain.
  ///   - `dgamma(mu)`  : free gluon leg γ^μ — `mu` is the OPEN Lorentz id (contracts the projector)
  ///   - `dslash(vlc)`: slashed/dressed propagator γ·p, `p = Σ coeff·fund(vid)` (mirrors @ref Elem::vlc)
  ///   - `dg5()`      : γ5 marker (chiral trace → Levi-Civita; staged separately)
  ///   - `dcomm(...)` : the BARE commutator `[A,B] = A·B − B·A` of two gammas as ONE chain token.
  ///                           Each leg A,B is independently a FREE leg (open Lorentz id) or a SLASH
  ///                           (momentum lin. comb.). This is the quark-gluon-vertex struct-7 tensor: it
  ///                           keeps the antisymmetric γ-pair folded so the commutator is never
  ///                           distributed into two separate traces. The σ^{μν}=(i/2)[γ^μ,γ^ν]
  ///                           normalization (the i/2 and any sign) lives in the emitted SCALAR, not the
  ///                           token — the front-end's Plus is already a bare bracket.
  struct DFac {
    enum Kind { Gamma, Gamma5, Slash, Comm, LoopSep };
    const Kind kind = Gamma;
    const int mu = -1; ///< kind 0: open Lorentz id; kind 3: leg-A FREE id (or -1 ⇒ leg-A is a slash, use `vlc`)
    const std::vector<std::pair<double, int>>
        vlc;           ///< kind 1: momentum lin. comb.; kind 3: leg-A slash momentum (when mu<0)
    const int nu = -1; ///< kind 3: leg-B FREE id (or -1 ⇒ leg-B is a slash, use `vlc2`)
    const std::vector<std::pair<double, int>> vlc2; ///< kind 3: leg-B slash momentum (when nu<0)
  };
  /// @brief A closed gamma chain in trace order (loop closes implicitly).
  using DiracNet = std::vector<DFac>;

  inline DFac dgamma(int muId) { return {DFac::Gamma, muId, {}, -1, {}}; } ///< free gluon leg γ^μ
  /// @brief Boundary between two INDEPENDENT closed spinor loops in one component (e.g. a quark loop
  ///        and the projection-closed external line, tied together only by gluon propagators). The
  ///        contraction traces each loop separately and multiplies the resulting Lorentz tensors,
  ///        contracting their shared gluon legs via the Lorentz net — NOT a Wick pairing across loops.
  inline DFac dloopsep() { return {DFac::LoopSep, -1, {}, -1, {}}; }
  inline DFac dslash(std::vector<std::pair<double, int>> vlc) { return {DFac::Slash, -1, std::move(vlc), -1, {}}; } ///< γ·p
  inline DFac dg5() { return {DFac::Gamma5, -1, {}, -1, {}}; }                                                      ///< γ5
  /// @brief bare commutator `[γ^μ, γ^ν] = γ^μγ^ν − γ^νγ^μ`, both legs FREE (open Lorentz ids μ,ν).
  inline DFac dcomm(int muId, int nuId) { return {DFac::Comm, muId, {}, nuId, {}}; }
  /// @brief bare commutator `[A̸, B̸]`, both legs SLASHED with momenta A,B (struct-7 external projector).
  inline DFac dcomm_ss(std::vector<std::pair<double, int>> a, std::vector<std::pair<double, int>> b)
  {
    return {DFac::Comm, -1, std::move(a), -1, std::move(b)};
  }
  /// @brief bare commutator `[γ^μ, B̸]`, leg-A FREE (gluon id μ), leg-B SLASHED with momentum B (loop vertex σ^{μν}B_ν).
  inline DFac dcomm_fs(int muId, std::vector<std::pair<double, int>> b)
  {
    return {DFac::Comm, muId, {}, -1, std::move(b)};
  }
  /// @brief bare commutator `[A̸, γ^ν]`, leg-A SLASHED with momentum A, leg-B FREE (gluon id ν).
  inline DFac dcomm_sf(std::vector<std::pair<double, int>> a, int nuId)
  {
    return {DFac::Comm, -1, std::move(a), nuId, {}};
  }

  namespace dirac_detail
  {

    /// @brief A vector leg carrying a full momentum linear combination `vlc` on Lorentz index `lbl`.
    inline NetVal vec_lc(int lbl, const std::vector<std::pair<double, int>> &vlc)
    {
      return {PTerm{Cx{1, 0}, {Elem{Elem::Vector, lbl, -1, -1, -1, vlc, 0, 0}}}};
    }

    /// @brief The pairing factor `g(a,b)` between two gamma tokens. Slash–slash uses a fresh shared
    ///        Lorentz label so the numeric contraction extracts the scalar product.
    inline NetVal pair_factor(const DFac &a, const DFac &b, int &fresh)
    {
      const bool af = (a.kind == DFac::Gamma), bf = (b.kind == DFac::Gamma);
      if (af && bf) return met(a.mu, b.mu);    // δ^{μa μb}
      if (af && !bf) return vec_lc(a.mu, b.vlc); // p_b^{μa}
      if (!af && bf) return vec_lc(b.mu, a.vlc); // p_a^{μb}
      const int sharedLbl = fresh++;             // slash–slash → p_a · p_b (sum over a fresh shared index)
      return contract(vec_lc(sharedLbl, a.vlc), vec_lc(sharedLbl, b.vlc));
    }

    /// @brief Wick pairing recursion for a token list of FREE/SLASH gammas (no γ5). Returns the trace
    ///        divided by 4 (the `tr(1)=4` is applied by @ref dirac_value).
    ///
    /// Pin the first token and pair it with each later token j: tr = Σ_{j≥1} (−1)^{j+1} g(t0,tj)·tr(rest),
    /// where `rest` is the chain with t0 and tj removed and g(·,·) is @ref pair_factor. Base case: the
    /// empty chain traces to 1 (= tr(1)/4).
    inline NetVal trace_rec(const std::vector<DFac> &tokens, int &fresh)
    {
      if (tokens.empty()) return {PTerm{Cx{1, 0}, {}}}; // empty trace: tr(1)/4 = 1 (tr(1)=4 applied by dirac_value)
      NetVal acc;                                       // empty == 0
      for (std::size_t j = 1; j < tokens.size(); ++j) {
        std::vector<DFac> rest;
        rest.reserve(tokens.size() - 2);
        for (std::size_t k = 1; k < tokens.size(); ++k)
          if (k != j) rest.push_back(tokens[k]);
        const double sgn = (j % 2 == 1) ? 1.0 : -1.0; // (−1)^{j+1}, 0-indexed (matches gammaTraceSum)
        NetVal pairFac = pair_factor(tokens[0], tokens[j], fresh);
        NetVal subTrace = trace_rec(rest, fresh);
        acc = add(std::move(acc), scale(sgn, contract(pairFac, subTrace)));
      }
      return acc;
    }

  } // namespace dirac_detail

  /// @brief Contract a closed gamma chain into a Lorentz @ref NetVal over its free legs.
  /// @param chain          the trace-ordered tokens (free legs / slashes; γ5 not yet handled here).
  /// @param firstFreeLabel a Lorentz id strictly above every label the surrounding component uses, so
  ///                       the fresh slash–slash labels never collide with the free legs / projector.
  /// @return the trace as a `NetVal` (empty == structural zero, e.g. an odd gamma count).
  inline NetVal dirac_value(const DiracNet &chain, int firstFreeLabel)
  {
    std::size_t ng = 0;
    for (const DFac &d : chain)
      if (d.kind != DFac::Gamma5) ++ng; // count gammas (free + slash)
    if (ng % 2 == 1) return {};         // odd chain → 0
    int fresh = firstFreeLabel;
    NetVal r = dirac_detail::trace_rec(chain, fresh); // stage 1: no γ5 tokens
    return scale(Cx{4.0, 0}, std::move(r));               // tr(1) = 4
  }

} // namespace numtracer::network
