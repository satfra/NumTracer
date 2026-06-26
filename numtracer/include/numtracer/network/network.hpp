/// @file network.hpp
/// @brief Lorentz network values: the symbolic Lorentz network (`NetVal`) and the
///        builders (`leaf`/`vec`/`met`/`proj`/`scale`/`add`/`contract`) that assemble it.
///
/// A network is a **sum of products** of metrics `δ_{μν}`, vectors `a^μ` (each carrying a
/// linear combination of momenta on one index — the eager-summation handle), and transverse
/// projectors `P(l)_{μν}`. The generator builds a diagram's Lorentz part with these builders,
/// then contracts it numerically to a scalar polynomial (@ref numtracer::numeric, in
/// `numeric/numeric_contract.hpp`).
#pragma once

#include "numtracer/core/cx.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace numtracer::network
{

  inline constexpr int kDim = 4; ///< Lorentz dimension (Euclidean here).

  /// @brief splitmix64 finaliser (local copy so this header has no heavy dependency; used by the
  ///        codegen env's open-addressed symbol index — constexpr-friendly).
  constexpr std::uint64_t splitmix64_finalise(std::uint64_t x)
  {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
  }
  /// @brief boost-style hash combine of two 64-bit hashes.
  constexpr std::uint64_t hash_combine(std::uint64_t a, std::uint64_t b)
  {
    return splitmix64_finalise(a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2)));
  }

  /// @brief A flattened network factor: kind 0 = metric(a,b), 1 = vector(a=index, vlc),
  ///        2 = transverse projector `P_T(l)`, 3 = Levi-Civita `ε_{a b c d}` (the γ5 trace's
  ///        totally-antisymmetric tensor), 4 = longitudinal projector `P_L(l)`.
  ///
  /// A `vector` factor carries a **linear combination** of momenta on one Lorentz index `a`
  /// (`vlc` = list of `(coeff, vid)`), not a single momentum. This is the eager-summation handle:
  /// a vertex sub-sum like `2·p_A − p_B` on a shared index is one compound vector leaf, so
  /// @ref contract never distributes it into separate product terms — the distribution that made
  /// the A4 reduction explode (21840 terms/net → ~81). The combination is expanded only at
  /// sp-extraction, after the (now far fewer) union-finds.
  ///
  /// An `epsilon` factor (kind 3) carries **four** Lorentz indices `a,b,c,d` (the `c`/`d` fields are
  /// used only by this kind). A closed Dirac trace has at most one γ5, hence at most one ε per
  /// product term — so the reduction never has to expand an ε·ε product (which would be a metric
  /// determinant); each ε's four indices simply contract with four momenta to a single antisymmetric
  /// invariant `eps(p_a,p_b,p_c,p_d)`. See @ref numtracer::numeric for the contraction.
  struct Elem {
    enum Kind { Metric, Vector, Epsilon, ProjT, ProjL, ProjE, ProjM };
    Kind kind = Metric;
    int a = 0, b = 0;                        ///< indices (kinds 0/2/3); kind 3 also uses c,d below
    int vid = -1;                            ///< proj's loop momentum `l` (kind 2)
    int inv = -1;                            ///< inverse env id (proj)
    std::vector<std::pair<double, int>> vlc; ///< vector linear combination `Σ coeff·vec(vid)` (kind 1)
    int c = 0, d = 0; ///< ε's 3rd/4th Lorentz indices (kind 3 only; default 0 leaves kinds 0/1/2 unchanged)
    int invS = -1;    ///< spatial inverse env id `1/|l⃗|²` (finite-T electric/magnetic projectors only)
  };

  // ---- network value: a sum of products (built by scale / add / contract) -----

  /// @brief One product term of a network: `coeff * prod(e)`.
  struct PTerm {
    Cx coeff{1, 0};
    std::vector<Elem> e;
  };
  /// @brief A network as a sum of product terms.
  using NetVal = std::vector<PTerm>;

  /// @brief A single-factor network (one product term, coefficient 1).
  inline constexpr NetVal leaf(Elem el) { return {PTerm{Cx{1, 0}, {el}}}; }
  /// @brief A vector leg `vid` on Lorentz index `Lbl`.
  inline constexpr NetVal vec(int Lbl, int Vid) { return leaf({Elem::Vector, Lbl, -1, -1, -1, {{1.0, Vid}}}); }
  /// @brief A metric `δ_{Mu Nu}`.
  inline constexpr NetVal met(int Mu, int Nu) { return leaf({Elem::Metric, Mu, Nu, -1, -1, {}}); }
  /// @brief A transverse projector `P_T(l)_{Mu Nu} = δ − l_Mu l_Nu/l²`, `l` = vector `Lvid`,
  ///        `1/l²` = env id `Inv`.
  inline constexpr NetVal projT(int Mu, int Nu, int Lvid, int Inv)
  {
    return leaf({Elem::ProjT, Mu, Nu, Lvid, Inv, {}});
  }
  /// @brief A longitudinal projector `P_L(l)_{Mu Nu} = l_Mu l_Nu/l²`, `l` = vector `Lvid`,
  ///        `1/l²` = env id `Inv`.
  inline constexpr NetVal projL(int Mu, int Nu, int Lvid, int Inv)
  {
    return leaf({Elem::ProjL, Mu, Nu, Lvid, Inv, {}});
  }
  /// @brief A finite-T **electric** (time-like-transverse) projector `P_E = P_T − P_M`, `l` = vector
  ///        `Lvid`, `1/l²` = env id `Inv`, `1/|l⃗|²` = env id `InvS`. Heat-bath direction is component 0.
  inline constexpr NetVal projE(int Mu, int Nu, int Lvid, int Inv, int InvS)
  {
    return leaf({Elem::ProjE, Mu, Nu, Lvid, Inv, {}, 0, 0, InvS});
  }
  /// @brief A finite-T **magnetic** (spatial-transverse) projector `P_M_{ij}=δ_{ij}−l_i l_j/|l⃗|²`
  ///        (i,j spatial; `P_M_{0ν}=P_M_{μ0}=0`), `l` = vector `Lvid`, `1/|l⃗|²` = env id `InvS`.
  inline constexpr NetVal projM(int Mu, int Nu, int Lvid, int InvS)
  {
    return leaf({Elem::ProjM, Mu, Nu, Lvid, -1, {}, 0, 0, InvS});
  }
  /// @brief A Levi-Civita tensor `ε_{Mu Nu Rho Sig}` (the γ5 trace's antisymmetric tensor).
  inline constexpr NetVal epsilon(int Mu, int Nu, int Rho, int Sig)
  {
    return leaf({Elem::Epsilon, Mu, Nu, -1, -1, {}, Rho, Sig});
  }

  /// @brief Multiply a network by a scalar.
  inline constexpr NetVal scale(Cx c, NetVal x)
  {
    for (PTerm &t : x)
      t.coeff = t.coeff * c;
    return x;
  }
  inline constexpr NetVal scale(double c, NetVal x) { return scale(Cx{c, 0}, std::move(x)); }

  /// @brief Whether `nv` is a pure sum of vectors all on the same Lorentz index — i.e. a momentum
  ///        linear combination that can collapse to one compound-vector leaf (eager summation).
  inline constexpr bool is_vecsum(const NetVal &nv, int &idx)
  {
    bool first = true;
    for (const PTerm &t : nv) {
      if (t.e.size() != 1 || t.e[0].kind != Elem::Vector) return false;
      if (first) {
        idx = t.e[0].a;
        first = false;
      } else if (t.e[0].a != idx)
        return false;
    }
    return !nv.empty();
  }

  /// @brief Sum of networks (concatenate their terms). **Eager summation:** when both summands are
  ///        vector sums on the same index (a vertex's momentum sub-combination), they collapse into a
  ///        single compound-vector leaf instead of two product terms, so @ref contract never
  ///        distributes the combination (the A4 explosion). Genuine structure sums (a vertex's sum of
  ///        metric×vector tensors) are *not* collapsible and concatenate as before.
  inline constexpr NetVal add(NetVal a, const NetVal &b)
  {
    int ia = 0, ib = 0;
    if (is_vecsum(a, ia) && is_vecsum(b, ib) && ia == ib) {
      Elem c;
      c.kind = Elem::Vector;
      c.a = ia;
      for (const PTerm &t : a)
        for (const auto &pr : t.e[0].vlc)
          c.vlc.push_back({pr.first * t.coeff.re, pr.second});
      for (const PTerm &t : b)
        for (const auto &pr : t.e[0].vlc)
          c.vlc.push_back({pr.first * t.coeff.re, pr.second});
      return {PTerm{Cx{1, 0}, {c}}};
    }
    a.insert(a.end(), b.begin(), b.end());
    return a;
  }
  template <class... R> constexpr NetVal add(NetVal a, const NetVal &b, const R &...r)
  {
    return add(add(std::move(a), b), r...);
  }

  /// @brief Tensor product of networks (Cartesian over their terms).
  inline constexpr NetVal contract(const NetVal &a, const NetVal &b)
  {
    NetVal r;
    r.reserve(a.size() * b.size());
    for (const PTerm &ta : a)
      for (const PTerm &tb : b) {
        PTerm p;
        p.coeff = ta.coeff * tb.coeff;
        p.e.reserve(ta.e.size() + tb.e.size());
        p.e.insert(p.e.end(), ta.e.begin(), ta.e.end());
        p.e.insert(p.e.end(), tb.e.begin(), tb.e.end());
        r.push_back(std::move(p));
      }
    return r;
  }
  template <class... R> constexpr NetVal contract(const NetVal &a, const NetVal &b, const R &...r)
  {
    return contract(contract(a, b), r...);
  }

} // namespace numtracer::network
