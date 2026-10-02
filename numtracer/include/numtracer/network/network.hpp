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
#include "numtracer/core/config.hpp" // NT_THROW (exception-optional guard for -fno-exceptions builds)

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace numtracer::inline network
{

  /// @brief A flattened network factor, tagged by @ref Elem::Kind:
  ///   - `Metric`  — δ_{a b}
  ///   - `Vector`  — a linear combination `Σ coeff·vec(vid)` of momenta on index `a` (held in `vlc`)
  ///   - `Epsilon` — Levi-Civita ε_{a b c d} (the γ5 trace's antisymmetric tensor; the only kind using c,d)
  ///   - `ProjT` / `ProjL` — transverse / longitudinal projectors `P_T(l)`, `P_L(l)`
  ///   - `ProjE` / `ProjM` — finite-T electric / magnetic projectors (use `invS`)
  ///
  /// A `Vector` factor carries a linear combination, not a single momentum, so a vertex sub-sum like
  /// `2·p_A − p_B` on a shared index stays one compound leaf and @ref contract never distributes it
  /// into separate terms — the distribution that made the A4 reduction explode (21840 terms/net → ~81).
  /// It is expanded only at scalar-product extraction, after the (now far fewer) union-finds.
  ///
  /// A closed Dirac trace has at most one γ5, hence at most one `Epsilon` per term, so the reduction
  /// never expands an ε·ε product; each ε's four indices contract with four momenta to a single
  /// invariant `eps(p_a,p_b,p_c,p_d)`. See @ref numtracer::numeric for the contraction.
  struct Elem {
    enum Kind { Metric, Vector, Epsilon, ProjT, ProjL, ProjE, ProjM };
    Kind kind = Metric;
    // Field docs name the VARIANT, never its enum ordinal: ordinals shift whenever `Kind` grows.
    // The ints sit together ahead of `vlc` so the struct packs into 56 B, not 64. Build it with
    // designated initializers only: a positional `Elem{…}` would silently mis-bind on a reorder.
    int a = 0, b = 0;                        ///< Lorentz index ids (Metric, Epsilon, and every projector)
    int vid = -1;                            ///< the projector's momentum `l` (ProjT / ProjL / ProjE / ProjM)
    int inv = -1;                            ///< inverse env id `1/l²` (projectors)
    int c = 0, d = 0; ///< ε's 3rd/4th Lorentz index ids (Epsilon only; the default 0 leaves every other kind unchanged)
    int invS = -1;    ///< spatial inverse env id `1/|l⃗|²` (finite-T electric/magnetic projectors only)
    std::vector<std::pair<double, int>> vlc{}; ///< vector linear combination `Σ coeff·vec(vid)` (Vector)
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
  inline NetVal leaf(Elem el) { return {PTerm{Cx{1, 0}, {el}}}; }
  /// @brief A vector leg `vid` on Lorentz index `Lbl`.
  inline NetVal vec(int Lbl, int Vid)
  {
    return leaf({.kind = Elem::Vector, .a = Lbl, .b = -1, .vid = -1, .inv = -1, .vlc = {{1.0, Vid}}});
  }
  /// @brief A metric `δ_{Mu Nu}`.
  inline NetVal met(int Mu, int Nu) { return leaf({.kind = Elem::Metric, .a = Mu, .b = Nu, .vid = -1, .inv = -1}); }
  /// @brief A transverse projector `P_T(l)_{Mu Nu} = δ − l_Mu l_Nu/l²`, `l` = vector `Lvid`,
  ///        `1/l²` = env id `Inv`.
  inline NetVal projT(int Mu, int Nu, int Lvid, int Inv)
  {
    return leaf({.kind = Elem::ProjT, .a = Mu, .b = Nu, .vid = Lvid, .inv = Inv});
  }
  /// @brief A longitudinal projector `P_L(l)_{Mu Nu} = l_Mu l_Nu/l²`, `l` = vector `Lvid`,
  ///        `1/l²` = env id `Inv`.
  inline NetVal projL(int Mu, int Nu, int Lvid, int Inv)
  {
    return leaf({.kind = Elem::ProjL, .a = Mu, .b = Nu, .vid = Lvid, .inv = Inv});
  }
  /// @brief A finite-T **electric** (time-like-transverse) projector `P_E = P_T − P_M`, `l` = vector
  ///        `Lvid`, `1/l²` = env id `Inv`, `1/|l⃗|²` = env id `InvS`. Heat-bath direction is component 0.
  inline NetVal projE(int Mu, int Nu, int Lvid, int Inv, int InvS)
  {
    return leaf({.kind = Elem::ProjE, .a = Mu, .b = Nu, .vid = Lvid, .inv = Inv, .invS = InvS});
  }
  /// @brief A finite-T **magnetic** (spatial-transverse) projector `P_M_{ij}=δ_{ij}−l_i l_j/|l⃗|²`
  ///        (i,j spatial; `P_M_{0ν}=P_M_{μ0}=0`), `l` = vector `Lvid`, `1/|l⃗|²` = env id `InvS`.
  inline NetVal projM(int Mu, int Nu, int Lvid, int InvS)
  {
    return leaf({.kind = Elem::ProjM, .a = Mu, .b = Nu, .vid = Lvid, .inv = -1, .invS = InvS});
  }
  /// @brief A Levi-Civita tensor `ε_{Mu Nu Rho Sig}` (the γ5 trace's antisymmetric tensor).
  inline NetVal epsilon(int Mu, int Nu, int Rho, int Sig)
  {
    return leaf({.kind = Elem::Epsilon, .a = Mu, .b = Nu, .vid = -1, .inv = -1, .c = Rho, .d = Sig});
  }

  /// @brief Multiply a network by a scalar.
  inline NetVal scale(Cx c, NetVal x)
  {
    for (PTerm &t : x)
      t.coeff = t.coeff * c;
    return x;
  }
  inline NetVal scale(double c, NetVal x) { return scale(Cx{c, 0}, std::move(x)); }

  /// @brief Whether `nv` is a pure sum of vectors all on the same Lorentz index — i.e. a momentum
  ///        linear combination that can collapse to one compound-vector leaf (eager summation).
  inline bool is_vecsum(const NetVal &nv, int &idx)
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
  ///        metric×vector tensors) are *not* collapsible and concatenate.
  inline NetVal add(NetVal a, const NetVal &b)
  {
    int ia = 0, ib = 0;
    if (is_vecsum(a, ia) && is_vecsum(b, ib) && ia == ib) {
      Elem c;
      c.kind = Elem::Vector;
      c.a = ia;
      // `vlc` coefficients are real by design (Elem::vlc is std::pair<double,int>): a momentum
      // linear combination has real weights. Fold each term's scalar coefficient into them — and
      // refuse a complex coefficient rather than silently dropping its imaginary part, which would
      // corrupt the contraction (the trap is a complex `scale` applied to a momentum vecsum
      // before this `add`).
      auto absorb = [&c](const NetVal &nv) {
        for (const PTerm &t : nv) {
          if (t.coeff.im != 0.0)
            NT_THROW(std::runtime_error, "network::add: complex coefficient on a vector-sum term "
                                         "(vlc weights are real-only)");
          for (const auto &pr : t.e[0].vlc)
            c.vlc.push_back({pr.first * t.coeff.re, pr.second});
        }
      };
      absorb(a);
      absorb(b);
      return {PTerm{Cx{1, 0}, {c}}};
    }
    a.insert(a.end(), b.begin(), b.end());
    return a;
  }
  template <class... R> NetVal add(NetVal a, const NetVal &b, const R &...r)
  {
    return add(add(std::move(a), b), r...);
  }

  /// @brief Tensor product of networks (Cartesian over their terms).
  inline NetVal contract(const NetVal &a, const NetVal &b)
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
  template <class... R> NetVal contract(const NetVal &a, const NetVal &b, const R &...r)
  {
    return contract(contract(a, b), r...);
  }

} // namespace numtracer::network
