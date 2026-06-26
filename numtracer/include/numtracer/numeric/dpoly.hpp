/// @file dpoly.hpp
/// @brief A **dressing polynomial**: a sum `Σ_m (dressing-monomial_m) · (kinematic MPoly_m)` carrying
///        runtime dressing/regulator calls as *symbolic atoms*, so a Feynman diagram whose propagator
///        numerators are dressed structure sums (e.g. `Mq·δ + Z(p)·γ·p`) collapses to **one** trace
///        instead of distributing into `2^D` separate diagrams (one per dressing combination).
///
/// This is the exact analogue of how @ref MPoly already carries `1/k²` inverse atoms in its monomial
/// key: a dressing monomial (@ref DMono) is a sorted multiset of *dressing-atom ids*, and a `DPoly`
/// maps each distinct dressing monomial to the (collected) kinematic @ref MPoly it multiplies. The
/// `MPoly` type itself is **unchanged** — `DPoly` is a thin wrapper whose arithmetic reuses
/// `MPoly::operator*`/`operator+` verbatim, so the numeric backend's hot path is untouched and only
/// diagrams that genuinely carry a dressed structure sum pay for the dressing layer.
///
/// At lowering (@ref numtracer::numeric::to_genprog) each dressing atom becomes a kind-2 `dress`
/// env leaf (@ref network::GlobalEnv::dr_id) — an opaque runtime value the kernel evaluates once,
/// exactly like an `inv(1/k²)` leaf — and the shared CSE/Horner pass collects the dressing factors
/// across monomials (FormTracer parity).
#pragma once

#include "numtracer/numeric/mpoly.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace numtracer::numeric
{

  /// @brief A dressing monomial: a sorted multiset of dressing-atom ids (the product of runtime
  ///        dressing/regulator calls multiplying one kinematic term). Empty == the un-dressed term.
  using DMono = std::vector<int>;

  /// @brief Canonicalise a dressing monomial (sort so equal products compare equal / collect).
  inline DMono dmono_sorted(DMono d)
  {
    std::sort(d.begin(), d.end());
    return d;
  }

  /// @brief Merge two sorted dressing monomials (multiset union) — the dressing analogue of the
  ///        atom-multiset merge in `MPoly::operator*` (`mpoly.hpp`).
  inline DMono dmono_merge(const DMono &a, const DMono &b)
  {
    DMono r;
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

  /// @brief A dressing polynomial: sorted-by-@ref DMono, like terms combined, no empty `MPoly` coeffs.
  struct DPoly {
    int nsym = 0;
    std::vector<std::pair<DMono, MPoly>> t; ///< sorted by DMono; each MPoly is non-empty

    DPoly() = default;
    explicit DPoly(int ns) : nsym(ns) {}

    bool empty() const { return t.empty(); }
    int size() const { return static_cast<int>(t.size()); }

    /// A `DPoly` that is just a single un-dressed kinematic polynomial (empty dressing monomial).
    /// The no-dressing case: lowering this is byte-for-byte the plain-`MPoly` path.
    static DPoly fromMPoly(const MPoly &p)
    {
      DPoly d(p.nsym);
      if (!p.empty()) d.t.push_back({DMono{}, p});
      return d;
    }

    /// Accumulate `p` into the coefficient of dressing monomial `d` (kept sorted; `d` already sorted).
    /// Drops the term if the resulting `MPoly` is empty (full cancellation).
    void add(const DMono &d, const MPoly &p)
    {
      if (p.empty()) return;
      auto it = std::lower_bound(t.begin(), t.end(), d,
                                 [](const std::pair<DMono, MPoly> &a, const DMono &k) { return a.first < k; });
      if (it != t.end() && it->first == d) {
        it->second = it->second + p;
        if (it->second.empty()) t.erase(it);
      } else {
        t.insert(it, {d, p});
      }
    }
  };

  inline DPoly operator+(const DPoly &a, const DPoly &b)
  {
    if (a.t.empty()) return b;
    if (b.t.empty()) return a;
    DPoly r(a.nsym ? a.nsym : b.nsym);
    r.t.reserve(a.t.size() + b.t.size());
    std::size_t i = 0, j = 0;
    while (i < a.t.size() && j < b.t.size()) {
      if (a.t[i].first < b.t[j].first)
        r.t.push_back(a.t[i++]);
      else if (b.t[j].first < a.t[i].first)
        r.t.push_back(b.t[j++]);
      else {
        MPoly s = a.t[i].second + b.t[j].second;
        if (!s.empty()) r.t.push_back({a.t[i].first, std::move(s)});
        ++i;
        ++j;
      }
    }
    while (i < a.t.size())
      r.t.push_back(a.t[i++]);
    while (j < b.t.size())
      r.t.push_back(b.t[j++]);
    return r;
  }

  /// Product: merge dressing monomials, multiply the kinematic `MPoly` coefficients (reusing
  /// `MPoly::operator*` verbatim), and collect. Provided for composability/testing; the contraction
  /// path builds a `DPoly` by accumulation (@ref DPoly::add) rather than multiplying two `DPoly`s.
  inline DPoly operator*(const DPoly &a, const DPoly &b)
  {
    const int ns = a.nsym ? a.nsym : b.nsym;
    DPoly r(ns);
    if (a.t.empty() || b.t.empty()) return r;
    for (const auto &[da, pa] : a.t)
      for (const auto &[db, pb] : b.t)
        r.add(dmono_merge(da, db), pa * pb);
    return r;
  }

  /// @brief Scale every kinematic coefficient by a complex constant (the group colour fold / sub-term
  ///        scalar). Empty terms drop out.
  inline DPoly scaleCx(const DPoly &a, Cx c)
  {
    DPoly r(a.nsym);
    if (c.re == 0 && c.im == 0) return r;
    const MPoly cc = MPoly::constant(a.nsym, c);
    r.t.reserve(a.t.size());
    for (const auto &[d, mp] : a.t) {
      MPoly s = mp * cc;
      if (!s.empty()) r.t.push_back({d, std::move(s)});
    }
    return r;
  }

  /// @brief Numeric evaluation (validation only): `x[i]` = user symbol i; `atomVal[aid]` = value of
  ///        `1/D_aid`; `drVal[id]` = value of dressing atom `id`. Equals the distributed sum
  ///        `Σ_combos (∏ dressings) · (kinematic value)`.
  inline Cx eval(const DPoly &p, const std::vector<double> &x, const std::vector<double> &atomVal,
                 const std::vector<double> &drVal)
  {
    Cx s{0, 0};
    for (const auto &[d, mp] : p.t) {
      double dm = 1.0;
      for (int id : d)
        dm *= drVal[id];
      const Cx kv = eval(mp, x, atomVal);
      s = s + Cx{kv.re * dm, kv.im * dm};
    }
    return s;
  }

} // namespace numtracer::numeric
