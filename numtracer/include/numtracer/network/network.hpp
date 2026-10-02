/// @file network.hpp
/// @brief The Lorentz network: index labels (@ref numtracer::LorentzIndex), momenta
///        (@ref numtracer::Momentum), and the network value (@ref numtracer::LorentzNet) with its
///        builders (`metric`, `vec`, `projT`/`projL`/`projE`/`projM`, `epsilon`) and algebra
///        (`*` = tensor product, `+` = sum, scalar multiples).
///
/// A network is a **sum of products** of metrics `δ_{μν}`, vectors `k_μ`, projectors `P(k)_{μν}` and
/// Levi-Civita tensors. Two factors are contracted exactly where they share a @ref LorentzIndex.
/// A @ref Frame contracts a finished network (together with a Dirac chain) to a scalar polynomial.
#pragma once

#include "numtracer/core/cx.hpp"
#include "numtracer/core/config.hpp" // NT_THROW (exception-optional guard for -fno-exceptions builds)

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace numtracer::inline network
{

  /// @brief A Lorentz index label. Two factors are summed over an index exactly when they carry the
  ///        same label (Einstein convention).
  ///
  /// Get fresh, distinct labels from @ref Frame::indices — `auto [mu, nu] = F.indices<2>();`. The
  /// constructor from a raw integer is `explicit` on purpose: it exists for generated code, which
  /// numbers its labels itself, and cannot be reached by accident (an integer, a vector id or an
  /// SU(N) label does not silently become a Lorentz index).
  struct LorentzIndex {
    int id;
    constexpr explicit LorentzIndex(int i) : id(i) {}
  };

  /// @brief A real linear combination of frame momenta, `Σ coeff·(momentum vid)`: the storage behind
  ///        @ref Momentum, and what the engine reads.
  using Vlc = std::vector<std::pair<double, int>>;

  /// @brief A momentum: a real linear combination of the momenta a @ref Frame declared.
  ///
  /// Obtained from @ref Frame::momentum and combined with `+`, `-` and real multiples, so `q = l - p`
  /// reads like the formula. The combination is kept sorted by momentum, with like terms merged and
  /// zero coefficients dropped.
  struct Momentum {
    Vlc lc; ///< `Σ coeff·(momentum vid)` as `{coeff, vid}` pairs
  };

  namespace mdetail
  {
    /// Sort by momentum id, merge equal ids, drop zero coefficients.
    inline Momentum canonical(Vlc lc)
    {
      std::sort(lc.begin(), lc.end(), [](const auto &x, const auto &y) { return x.second < y.second; });
      Vlc out;
      for (const auto &[c, v] : lc) {
        if (!out.empty() && out.back().second == v)
          out.back().first += c;
        else
          out.push_back({c, v});
      }
      out.erase(std::remove_if(out.begin(), out.end(), [](const auto &t) { return t.first == 0.0; }), out.end());
      return {std::move(out)};
    }
  } // namespace mdetail

  inline Momentum operator+(Momentum a, const Momentum &b)
  {
    a.lc.insert(a.lc.end(), b.lc.begin(), b.lc.end());
    return mdetail::canonical(std::move(a.lc));
  }
  inline Momentum operator*(double c, Momentum a)
  {
    for (auto &t : a.lc) t.first *= c;
    return mdetail::canonical(std::move(a.lc));
  }
  inline Momentum operator*(Momentum a, double c) { return c * std::move(a); }
  inline Momentum operator-(Momentum a) { return -1.0 * std::move(a); }
  inline Momentum operator-(Momentum a, const Momentum &b) { return std::move(a) + (-b); }

  /// @brief One factor of a Lorentz network, tagged by @ref LorentzFactor::Kind:
  ///   - `Metric`  — δ_{a b}
  ///   - `Vector`  — the momentum `vlc` on index `a`
  ///   - `Epsilon` — Levi-Civita ε_{a b c d} (the only kind using c, d)
  ///   - `ProjT` / `ProjL` — transverse / longitudinal projectors `P_T(k)_{ab}`, `P_L(k)_{ab}`
  ///   - `ProjE` / `ProjM` — finite-T electric / magnetic projectors (these use `atomS`)
  ///
  /// You normally build factors with the builders below, not by hand. This is plain data; the
  /// generated code fills it with designated initializers.
  ///
  /// A projector's momentum `k` is `vid` when it is a single frame momentum (`vlc` empty — the common
  /// case, which keeps the factor allocation-free), else the linear combination `vlc`. A vector always
  /// carries `vlc`, so a vertex sub-sum like `2·p_A − p_B` on one index stays one compound factor and
  /// @ref mul never distributes it into separate terms.
  struct LorentzFactor {
    enum Kind { Metric, Vector, Epsilon, ProjT, ProjL, ProjE, ProjM };
    Kind kind = Metric;
    // Field docs name the VARIANT, never its enum ordinal: ordinals shift whenever `Kind` grows.
    // The ints sit together ahead of `vlc` so the struct packs into 56 B, not 64. Build it with
    // designated initializers only: a positional `LorentzFactor{…}` would silently mis-bind on a reorder.
    int a = 0, b = 0; ///< Lorentz index ids (Metric, Epsilon, and every projector)
    int vid = -1;     ///< a projector's momentum when it is the single frame momentum `vid` (`vlc` empty)
    int atom = -1;    ///< id of the projector's `1/k²` (ProjT / ProjL / ProjE); -1 = let the Frame assign it
    int c = 0, d = 0; ///< ε's 3rd/4th Lorentz index ids (Epsilon only)
    int atomS = -1;   ///< id of the spatial `1/|k⃗|²` (ProjE / ProjM); -1 = let the Frame assign it
    Vlc vlc{};        ///< a vector's momentum, or a projector's when it is a linear combination

    /// Whether this is one of the four projector kinds.
    bool is_projector() const { return kind == ProjT || kind == ProjL || kind == ProjE || kind == ProjM; }
  };
  static_assert(sizeof(LorentzFactor) <= 56, "LorentzFactor grew: it is copied per term on the hot path");

  // ---- network value: a sum of products -----

  /// @brief One product term of a network: `coeff * prod(e)`.
  struct LorentzTerm {
    Cx coeff{1, 0};
    std::vector<LorentzFactor> e;
  };
  /// @brief A Lorentz network: a sum of product terms. An empty network is the scalar 1 when handed
  ///        to a contraction.
  using LorentzNet = std::vector<LorentzTerm>;

  /// @brief A single-factor network (one product term, coefficient 1).
  inline LorentzNet leaf(LorentzFactor el) { return {LorentzTerm{Cx{1, 0}, {std::move(el)}}}; }

  namespace mdetail
  {
    /// A projector factor on momentum @p k: a single frame momentum rides `vid`, anything else `vlc`.
    inline LorentzFactor projector(LorentzFactor::Kind kind, LorentzIndex mu, LorentzIndex nu, const Momentum &k,
                                   int atom, int atomS)
    {
      if (k.lc.empty())
        NT_THROW(std::invalid_argument, "projector on a zero momentum: its 1/k^2 is undefined");
      LorentzFactor f{.kind = kind, .a = mu.id, .b = nu.id, .atom = atom, .atomS = atomS};
      if (k.lc.size() == 1 && k.lc[0].first == 1.0)
        f.vid = k.lc[0].second;
      else
        f.vlc = k.lc;
      return f;
    }
  } // namespace mdetail

  /// @brief The metric `δ_{μν}`.
  inline LorentzNet metric(LorentzIndex mu, LorentzIndex nu)
  {
    return leaf({.kind = LorentzFactor::Metric, .a = mu.id, .b = nu.id});
  }
  /// @brief The vector `k_μ`.
  inline LorentzNet vec(LorentzIndex mu, const Momentum &k)
  {
    return leaf({.kind = LorentzFactor::Vector, .a = mu.id, .b = -1, .vlc = k.lc});
  }
  /// @brief The transverse projector `P_T(k)_{μν} = δ_{μν} − k_μ k_ν/k²`.
  /// @param atom id of its `1/k²`; leave it at -1 and the @ref Frame assigns one.
  inline LorentzNet projT(LorentzIndex mu, LorentzIndex nu, const Momentum &k, int atom = -1)
  {
    return leaf(mdetail::projector(LorentzFactor::ProjT, mu, nu, k, atom, -1));
  }
  /// @brief The longitudinal projector `P_L(k)_{μν} = k_μ k_ν/k²`.
  inline LorentzNet projL(LorentzIndex mu, LorentzIndex nu, const Momentum &k, int atom = -1)
  {
    return leaf(mdetail::projector(LorentzFactor::ProjL, mu, nu, k, atom, -1));
  }
  /// @brief The finite-T **electric** projector `P_E = P_T − P_M` (heat-bath direction = component 0).
  inline LorentzNet projE(LorentzIndex mu, LorentzIndex nu, const Momentum &k, int atom = -1, int atomS = -1)
  {
    return leaf(mdetail::projector(LorentzFactor::ProjE, mu, nu, k, atom, atomS));
  }
  /// @brief The finite-T **magnetic** projector `P_M_{ij} = δ_{ij} − k_i k_j/|k⃗|²` on the spatial
  ///        components (row and column 0 vanish).
  inline LorentzNet projM(LorentzIndex mu, LorentzIndex nu, const Momentum &k, int atomS = -1)
  {
    return leaf(mdetail::projector(LorentzFactor::ProjM, mu, nu, k, -1, atomS));
  }
  /// @brief The Levi-Civita tensor `ε_{μνρσ}`, `ε_{0123} = +1`.
  inline LorentzNet epsilon(LorentzIndex mu, LorentzIndex nu, LorentzIndex rho, LorentzIndex sigma)
  {
    return leaf({.kind = LorentzFactor::Epsilon, .a = mu.id, .b = nu.id, .c = rho.id, .d = sigma.id});
  }

  /// @brief Multiply a network by a scalar.
  inline LorentzNet scale(Cx c, LorentzNet x)
  {
    for (LorentzTerm &t : x)
      t.coeff = t.coeff * c;
    return x;
  }
  inline LorentzNet scale(double c, LorentzNet x) { return scale(Cx{c, 0}, std::move(x)); }

  /// @brief Whether `nv` is a pure sum of vectors all on the same Lorentz index — i.e. a momentum
  ///        linear combination that can collapse to one compound-vector leaf (eager summation).
  inline bool is_vecsum(const LorentzNet &nv, int &idx)
  {
    bool first = true;
    for (const LorentzTerm &t : nv) {
      if (t.e.size() != 1 || t.e[0].kind != LorentzFactor::Vector) return false;
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
  ///        single compound-vector leaf instead of two product terms, so @ref mul never
  ///        distributes the combination (the A4 explosion). Genuine structure sums (a vertex's sum of
  ///        metric×vector tensors) are *not* collapsible and concatenate.
  inline LorentzNet add(LorentzNet a, const LorentzNet &b)
  {
    int ia = 0, ib = 0;
    if (is_vecsum(a, ia) && is_vecsum(b, ib) && ia == ib) {
      LorentzFactor c;
      c.kind = LorentzFactor::Vector;
      c.a = ia;
      // `vlc` coefficients are real by design (Vlc holds real weights): a momentum
      // linear combination has real weights. Fold each term's scalar coefficient into them — and
      // refuse a complex coefficient rather than silently dropping its imaginary part, which would
      // corrupt the contraction (the trap is a complex `scale` applied to a momentum vecsum
      // before this `add`).
      auto absorb = [&c](const LorentzNet &nv) {
        for (const LorentzTerm &t : nv) {
          if (t.coeff.im != 0.0)
            NT_THROW(std::runtime_error, "add: complex coefficient on a vector-sum term "
                                         "(vlc weights are real-only)");
          for (const auto &pr : t.e[0].vlc)
            c.vlc.push_back({pr.first * t.coeff.re, pr.second});
        }
      };
      absorb(a);
      absorb(b);
      return {LorentzTerm{Cx{1, 0}, {c}}};
    }
    a.insert(a.end(), b.begin(), b.end());
    return a;
  }
  template <class... R> LorentzNet add(LorentzNet a, const LorentzNet &b, const R &...r)
  {
    return add(add(std::move(a), b), r...);
  }

  /// @brief Tensor product of networks (Cartesian over their terms): factors sharing a label are
  ///        contracted when the result is evaluated. Spelled `a * b` in user code.
  inline LorentzNet mul(const LorentzNet &a, const LorentzNet &b)
  {
    LorentzNet r;
    r.reserve(a.size() * b.size());
    for (const LorentzTerm &ta : a)
      for (const LorentzTerm &tb : b) {
        LorentzTerm p;
        p.coeff = ta.coeff * tb.coeff;
        p.e.reserve(ta.e.size() + tb.e.size());
        p.e.insert(p.e.end(), ta.e.begin(), ta.e.end());
        p.e.insert(p.e.end(), tb.e.begin(), tb.e.end());
        r.push_back(std::move(p));
      }
    return r;
  }
  template <class... R> LorentzNet mul(const LorentzNet &a, const LorentzNet &b, const R &...r)
  {
    return mul(mul(a, b), r...);
  }

  /// @brief `a * b`: the tensor product (@ref mul).
  inline LorentzNet operator*(const LorentzNet &a, const LorentzNet &b) { return mul(a, b); }
  /// @brief `a + b`: the sum (@ref add).
  inline LorentzNet operator+(LorentzNet a, const LorentzNet &b) { return add(std::move(a), b); }
  /// @brief A scalar multiple of a network.
  inline LorentzNet operator*(Cx c, LorentzNet x) { return scale(c, std::move(x)); }
  inline LorentzNet operator*(double c, LorentzNet x) { return scale(c, std::move(x)); }

} // namespace numtracer::network
