/// @file dressed_poly.hpp
/// @brief A **dressing polynomial**: a sum `Σ_m (dressing-monomial_m) · (kinematic Poly_m)` carrying
///        runtime dressing/regulator calls as *symbolic atoms*, so a Feynman diagram whose propagator
///        numerators are dressed structure sums (e.g. `Mq·δ + Z(p)·γ·p`) collapses to **one** trace
///        instead of distributing into `2^D` separate diagrams (one per dressing combination).
///
/// This is the exact analogue of how @ref Poly already carries `1/k²` inverse atoms in its monomial
/// key: a dressing monomial (@ref DressingMono) is a sorted multiset of *dressing-atom ids*, and a `DressedPoly`
/// maps each distinct dressing monomial to the (collected) kinematic @ref Poly it multiplies. The
/// `Poly` type itself is **unchanged** — `DressedPoly` is a thin wrapper whose arithmetic reuses
/// `Poly::operator*`/`operator+` verbatim, so the numeric backend's hot path is untouched and only
/// diagrams that genuinely carry a dressed structure sum pay for the dressing layer.
///
/// At lowering (@ref numtracer::network::to_genprog) each dressing atom becomes a `SymKind::dress`
/// env leaf (@ref network::GlobalEnv::dr_id) — an opaque runtime value the kernel evaluates once,
/// exactly like an `inv(1/k²)` leaf — and the shared CSE/Horner pass collects the dressing factors
/// across monomials (FormTracer parity).
#pragma once

#include "numtracer/poly/poly.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace numtracer::inline poly
{

  /// @brief A dressing monomial: a sorted multiset of dressing-atom ids (the product of runtime
  ///        dressing/regulator calls multiplying one kinematic term). Empty == the un-dressed term.
  using DressingMono = std::vector<int>;

  /// @brief Canonicalise a dressing monomial (sort so equal products compare equal / collect).
  inline DressingMono dressing_mono_sorted(DressingMono d)
  {
    std::sort(d.begin(), d.end());
    return d;
  }

  /// @brief Merge two sorted dressing monomials (multiset union) — the dressing analogue of the
  ///        atom-multiset merge in `Poly::operator*` (`poly.hpp`).
  inline DressingMono dressing_mono_merge(const DressingMono &a, const DressingMono &b)
  {
    DressingMono r;
    r.reserve(a.size() + b.size());
    std::size_t i = 0, j = 0;
    while (i < a.size() && j < b.size())
      r.push_back(a[i] <= b[j] ? a[i++] : b[j++]);
    while (i < a.size())
      r.push_back(a[i++]);
    while (j < b.size())
      r.push_back(b[j++]);
    return r;
  }

  // The `nsym`-carrying construction API is closed behind these friends, exactly as for @ref Poly:
  // @ref Frame is the sole user-facing path and @ref DressedPolyFactory the internal attorney for the
  // trusted contraction/generator-phase code. (`poly.hpp`, included above, already forward-declares
  // @ref Frame; declare the DressedPoly attorney here.)
  struct DressedPolyFactory;

  /// @brief A dressing polynomial: sorted-by-@ref DressingMono, like terms combined, no empty `Poly` coeffs.
  struct DressedPoly {
    int nsym = 0;
    std::vector<std::pair<DressingMono, Poly>> terms; ///< sorted by DressingMono; each Poly is non-empty

    // Sanctioned construction paths (see @ref Poly). The in-header DressedPoly arithmetic constructs its
    // results directly.
    friend class engine::Frame;
    friend struct DressedPolyFactory;
    friend DressedPoly operator+(const DressedPoly &a, const DressedPoly &b);
    friend DressedPoly operator*(const DressedPoly &a, const DressedPoly &b);
    friend DressedPoly scaleCx(const DressedPoly &a, Cx c);

    DressedPoly() = default;

  private:
    // Bare-`nsym` construction — reachable only through @ref Frame / @ref DressedPolyFactory (friends);
    // see @ref Poly. The empty default ctor above stays public.
    explicit DressedPoly(int ns) : nsym(ns) {}

    /// A `DressedPoly` that is just a single un-dressed kinematic polynomial (empty dressing monomial).
    /// The no-dressing case: lowering this is byte-for-byte the plain-`Poly` path.
    static DressedPoly fromPoly(const Poly &p)
    {
      DressedPoly d(p.nsym);
      if (!p.empty()) d.terms.push_back({DressingMono{}, p});
      return d;
    }

  public:
    bool empty() const { return terms.empty(); }
    int size() const { return static_cast<int>(terms.size()); }

    /// Accumulate `p` into the coefficient of dressing monomial `d` (kept sorted; `d` already sorted).
    /// Drops the term if the resulting `Poly` is empty (full cancellation).
    void add(const DressingMono &d, const Poly &p)
    {
      if (p.empty()) return;
      auto it =
          std::lower_bound(terms.begin(), terms.end(), d,
                           [](const std::pair<DressingMono, Poly> &a, const DressingMono &k) { return a.first < k; });
      if (it != terms.end() && it->first == d) {
        it->second = it->second + p;
        if (it->second.empty()) terms.erase(it);
      } else {
        terms.insert(it, {d, p});
      }
    }
  };

  /// @brief Internal attorney re-exposing the private @ref DressedPoly factories to the trusted engine code
  ///        (the dressed contraction / generator phases), mirroring @ref PolyFactory. Not public API.
  struct DressedPolyFactory {
    static DressedPoly zero(int ns) { return DressedPoly(ns); }
    static DressedPoly fromPoly(const Poly &p) { return DressedPoly::fromPoly(p); }
  };

  inline DressedPoly operator+(const DressedPoly &a, const DressedPoly &b)
  {
    if (a.terms.empty()) return b;
    if (b.terms.empty()) return a;
    DressedPoly r(a.nsym ? a.nsym : b.nsym);
    r.terms.reserve(a.terms.size() + b.terms.size());
    std::size_t i = 0, j = 0;
    while (i < a.terms.size() && j < b.terms.size()) {
      if (a.terms[i].first < b.terms[j].first)
        r.terms.push_back(a.terms[i++]);
      else if (b.terms[j].first < a.terms[i].first)
        r.terms.push_back(b.terms[j++]);
      else {
        Poly s = a.terms[i].second + b.terms[j].second;
        if (!s.empty()) r.terms.push_back({a.terms[i].first, std::move(s)});
        ++i;
        ++j;
      }
    }
    while (i < a.terms.size())
      r.terms.push_back(a.terms[i++]);
    while (j < b.terms.size())
      r.terms.push_back(b.terms[j++]);
    return r;
  }

  /// Product: merge dressing monomials, multiply the kinematic `Poly` coefficients (reusing
  /// `Poly::operator*` verbatim), and collect. Provided for composability/testing; the contraction
  /// path builds a `DressedPoly` by accumulation (@ref DressedPoly::add) rather than multiplying two `DressedPoly`s.
  inline DressedPoly operator*(const DressedPoly &a, const DressedPoly &b)
  {
    const int ns = a.nsym ? a.nsym : b.nsym;
    DressedPoly r(ns);
    if (a.terms.empty() || b.terms.empty()) return r;
    for (const auto &[da, pa] : a.terms)
      for (const auto &[db, pb] : b.terms)
        r.add(dressing_mono_merge(da, db), pa * pb);
    return r;
  }

  /// @brief Scale every kinematic coefficient by a complex constant (the group colour fold / sub-term
  ///        scalar). Empty terms drop out.
  inline DressedPoly scaleCx(const DressedPoly &a, Cx c)
  {
    DressedPoly r(a.nsym);
    if (c.re == 0 && c.im == 0) return r;
    r.terms.reserve(a.terms.size());
    for (const auto &[d, mp] : a.terms) {
      // Direct coefficient scaling instead of `mp * constant(c)`: bit-identical (see Poly::scaled),
      // without the scratch and sort. The emptiness guard is defensive: `mp` is non-empty and `c != 0`.
      Poly s = PolyFactory::scaled(mp, c);
      if (!s.empty()) r.terms.push_back({d, std::move(s)});
    }
    return r;
  }

  namespace pdetail
  {
    /// @brief Raw numeric evaluation of a dressed polynomial: `x[i]` = user symbol i; `atomVal[aid]` =
    ///        value of `1/D_aid`; `drVal[id]` = value of dressing atom `id`. Equals the distributed sum
    ///        `Σ_combos (∏ dressings) · (kinematic value)`. Prefer @ref Frame::eval.
    inline Cx eval(const DressedPoly &p, const std::vector<double> &x, const std::vector<double> &atomVal,
                   const std::vector<double> &drVal)
    {
      Cx s{0, 0};
      for (const auto &[d, mp] : p.terms) {
        double dm = 1.0;
        for (int id : d) {
          if (id < 0 || static_cast<std::size_t>(id) >= drVal.size())
            NT_THROW(std::invalid_argument, ("eval: the polynomial carries dressing " + std::to_string(id) +
                                             " but only " + std::to_string(drVal.size()) + " dressing values were given")
                                                .c_str());
          dm *= drVal[id];
        }
        const Cx kv = eval(mp, x, atomVal);
        s = s + Cx{kv.re * dm, kv.im * dm};
      }
      return s;
    }
  } // namespace pdetail

} // namespace numtracer::inline poly
