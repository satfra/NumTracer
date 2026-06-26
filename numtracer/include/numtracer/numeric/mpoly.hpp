/// @file mpoly.hpp
/// @brief Multivariate polynomial over an arbitrary user-symbol set, carrying surviving inverse
///        (`1/k²`) atoms in the monomial key. The arithmetic core of the **numeric contraction
///        backend** (task #22): contract a diagram as 4×4 spinor matrix products whose entries are
///        these polynomials, so γ/metric/projector stay numeric and only the user's symbolic
///        momentum data enters the result — bounded by the matrix structure, not the `(2n−1)!!`
///        sp-monomial blowup of the symbolic `reduce_product` path.
///
/// A `MPoly` is a map from monomial → complex coefficient. A monomial is
///   (exponent vector over the `nsym` user symbols, a sorted multiset of inverse-atom ids).
/// The atom ids ride along so a transverse projector's `INV(k) = 1/k²` factor is tracked exactly:
/// when a numerator monomial is divisible by an atom whose **denominator is a single monomial**
/// (e.g. a unit-direction loop `k²=l1²`), @ref divThroughMonomialAtoms cancels it; an atom whose
/// denominator is a genuine polynomial (a shifted line `k=l−q`, `k²` non-monomial) survives and is
/// later lowered to an `inv` env slot. This is the general bare-loop cancellation, not a special case.
///
/// Map-based (sizes are tens of monomials with frame inputs — the regime where this backend wins);
/// the symbol *meaning* lives only in the kernel's `fill`, so the engine is frame-agnostic.
///
/// This header is the polynomial type alone; the 4×4 spinor matrix of `MPoly` (`Mat4`) and the
/// γ/slash builders that consume it live in `numeric/spinor_mat.hpp`.
#pragma once

#include "numtracer/core/cx.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

namespace numtracer::numeric
{

  using numtracer::Cx;

  /// @brief A monomial: exponents over the user symbols plus a sorted multiset of inverse-atom ids.
  struct Mono {
    std::vector<int> e;     ///< length nsym, exponent of each user symbol
    std::vector<int> atoms; ///< sorted (with multiplicity) ids of surviving `1/k²` atoms
    bool operator<(const Mono &o) const
    {
      if (e != o.e) return e < o.e;
      return atoms < o.atoms;
    }
    bool operator==(const Mono &o) const { return e == o.e && atoms == o.atoms; }
  };

  /// @brief Multivariate polynomial: monomial → complex coefficient, over `nsym` user symbols.
  ///
  /// Storage is a **sorted, like-terms-combined `std::vector`** (not a `std::map`): the per-insert
  /// red-black-tree overhead of a map dominates for the large (10⁴+ monomial) pure-gauge polynomials,
  /// so `operator*` collects all `n·m` products into a scratch vector and **sort-collects once**
  /// (O(nm log nm), no per-term tree churn) and `operator+` is a linear merge. `t` stays sorted by
  /// @ref Mono and carries no zero coefficients, so iteration order is deterministic (reproducible
  /// kernel) and equality/lookup are binary searches.
  struct MPoly {
    int nsym = 0;
    std::vector<std::pair<Mono, Cx>> t; ///< sorted by Mono, like terms combined, no zeros

    MPoly() = default;
    explicit MPoly(int ns) : nsym(ns) {}

    /// Build from an unsorted scratch list of (monomial, coeff): sort then combine adjacent equals.
    static MPoly from_scratch(int ns, std::vector<std::pair<Mono, Cx>> s)
    {
      MPoly p(ns);
      std::sort(s.begin(), s.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
      p.t.reserve(s.size());
      for (auto &kv : s) {
        if (kv.second.re == 0 && kv.second.im == 0) continue;
        if (!p.t.empty() && p.t.back().first == kv.first) {
          p.t.back().second = p.t.back().second + kv.second;
          if (p.t.back().second.re == 0 && p.t.back().second.im == 0) p.t.pop_back();
        } else {
          p.t.push_back(std::move(kv));
        }
      }
      return p;
    }

    static MPoly constant(int ns, Cx c)
    {
      MPoly p(ns);
      if (!(c.re == 0 && c.im == 0)) p.t.push_back({Mono{std::vector<int>(ns, 0), {}}, c});
      return p;
    }
    /// The i-th user symbol (coefficient 1).
    static MPoly var(int ns, int i)
    {
      MPoly p(ns);
      std::vector<int> e(ns, 0);
      e[i] = 1;
      p.t.push_back({Mono{std::move(e), {}}, Cx{1, 0}});
      return p;
    }
    /// A single monomial `c · ∏ x_k^{e_k}` (no inverse atoms). Used by the generated component table.
    static MPoly mono(int ns, std::vector<int> e, Cx c)
    {
      MPoly p(ns);
      if (!(c.re == 0 && c.im == 0)) p.t.push_back({Mono{std::move(e), {}}, c});
      return p;
    }
    /// A bare inverse atom `1/D` (atom id `aid`), coefficient 1.
    static MPoly atom(int ns, int aid)
    {
      MPoly p(ns);
      p.t.push_back({Mono{std::vector<int>(ns, 0), {aid}}, Cx{1, 0}});
      return p;
    }

    int size() const { return t.size(); }
    bool empty() const { return t.empty(); }

    /// Insert/accumulate one term (keeps `t` sorted). O(n) shift — used only on the incremental paths
    /// (component builders); the hot `operator*`/`operator+`/reductions build whole vectors at once.
    void addTerm(const Mono &m, Cx c)
    {
      if (c.re == 0 && c.im == 0) return;
      auto it = std::lower_bound(t.begin(), t.end(), m,
                                 [](const std::pair<Mono, Cx> &a, const Mono &k) { return a.first < k; });
      if (it != t.end() && it->first == m) {
        it->second = it->second + c;
        if (it->second.re == 0 && it->second.im == 0) t.erase(it);
      } else {
        t.insert(it, {m, c});
      }
    }
  };

  inline MPoly operator+(const MPoly &a, const MPoly &b)
  {
    if (a.t.empty()) return b;
    if (b.t.empty()) return a;
    // Both operands carry terms, so their symbol spaces must agree: the merge below walks the two
    // exponent vectors slot-for-slot, and a mismatch would read out of bounds (a default-constructed
    // zero with nsym == 0 is allowed — it is caught by the empty checks above). Debug-only; compiles
    // out under NDEBUG.
    assert(a.nsym == b.nsym);
    const int ns = a.nsym ? a.nsym : b.nsym;
    MPoly r(ns);
    r.t.reserve(a.t.size() + b.t.size());
    std::size_t i = 0, j = 0;
    while (i < a.t.size() && j < b.t.size()) {
      if (a.t[i].first < b.t[j].first)
        r.t.push_back(a.t[i++]);
      else if (b.t[j].first < a.t[i].first)
        r.t.push_back(b.t[j++]);
      else {
        Cx s = a.t[i].second + b.t[j].second;
        if (!(s.re == 0 && s.im == 0)) r.t.push_back({a.t[i].first, s});
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
  inline MPoly operator-(const MPoly &a, const MPoly &b)
  {
    assert(a.t.empty() || b.t.empty() || a.nsym == b.nsym); // see operator+; debug-only
    MPoly nb(b.nsym);
    nb.t.reserve(b.t.size());
    for (const auto &kv : b.t)
      nb.t.push_back({kv.first, Cx{-kv.second.re, -kv.second.im}});
    return a + nb;
  }
  inline MPoly operator*(const MPoly &a, const MPoly &b)
  {
    const int ns = a.nsym ? a.nsym : b.nsym;
    if (a.t.empty() || b.t.empty()) return MPoly(ns);
    assert(a.nsym == b.nsym); // both carry terms ⇒ symbol spaces must match; see operator+ (debug-only)
    std::vector<std::pair<Mono, Cx>> s;
    s.reserve(a.t.size() * b.t.size());
    std::vector<int> e(ns);
    for (const auto &[ma, ca] : a.t)
      for (const auto &[mb, cb] : b.t) {
        for (int k = 0; k < ns; ++k)
          e[k] = ma.e[k] + mb.e[k];
        std::vector<int> at;
        at.reserve(ma.atoms.size() + mb.atoms.size());
        std::size_t i = 0, j = 0; // merge the two sorted atom multisets
        while (i < ma.atoms.size() && j < mb.atoms.size())
          at.push_back(ma.atoms[i] <= mb.atoms[j] ? ma.atoms[i++] : mb.atoms[j++]);
        while (i < ma.atoms.size())
          at.push_back(ma.atoms[i++]);
        while (j < mb.atoms.size())
          at.push_back(mb.atoms[j++]);
        s.push_back({Mono{e, std::move(at)}, ca * cb});
      }
    return MPoly::from_scratch(ns, std::move(s));
  }

  /// @brief Cancel each numerator monomial against any atom whose denominator is a single monomial.
  ///
  /// `atomDen[aid]` is the atom's denominator polynomial `k²`. If it is a single monomial
  /// `c·∏ x^d`, then one power of the atom `1/D` cancels whenever the numerator exponents dominate
  /// `d` componentwise: subtract `d`, divide the coefficient by `c`, drop one atom instance. Repeat
  /// per atom occurrence. Atoms with a non-monomial denominator (shifted lines) are left intact and
  /// survive into the lowering as `inv` env slots. Value-preserving and frame-agnostic.
  inline MPoly divThroughMonomialAtoms(const MPoly &p, const std::vector<MPoly> &atomDen)
  {
    std::vector<std::pair<Mono, Cx>> out;
    out.reserve(p.t.size());
    for (const auto &[m, c] : p.t) {
      std::vector<int> e = m.e;
      Cx coeff = c;
      std::vector<int> keep;
      keep.reserve(m.atoms.size());
      for (const auto &aid : m.atoms) {
        bool cancelled = false;
        if (aid >= 0 && aid < (int)atomDen.size()) {
          const MPoly &D = atomDen[aid];
          if (D.t.size() == 1) { // monomial denominator → may cancel
            const auto &dm = *D.t.begin();
            const std::vector<int> &d = dm.first.e;
            const Cx dc = dm.second;
            bool dominates = dm.first.atoms.empty();
            for (std::size_t k = 0; dominates && k < e.size(); ++k)
              if (e[k] < d[k]) dominates = false;
            if (dominates) {
              for (std::size_t k = 0; k < e.size(); ++k)
                e[k] -= d[k];
              // coeff /= dc   (complex division: z / w = z·conj(w) / |w|²)
              const double den = dc.re * dc.re + dc.im * dc.im;
              coeff = Cx{(coeff.re * dc.re + coeff.im * dc.im) / den, (coeff.im * dc.re - coeff.re * dc.im) / den};
              cancelled = true;
            }
          }
        }
        if (!cancelled) keep.push_back(aid);
      }
      // keep stays sorted because m.atoms was sorted and we only dropped elements
      if (!(coeff.re == 0 && coeff.im == 0)) out.push_back({Mono{std::move(e), std::move(keep)}, coeff});
    }
    return MPoly::from_scratch(p.nsym, std::move(out));
  }

  /// @brief Apply unit-vector constraints `Σ_μ Uμ² = 1` to reduce each unit group's LAST component to
  ///        power ≤ 1: `U_last² → 1 − Σ_{μ<last} Uμ²`. A `group` is the list of symbol indices of one
  ///        unit direction's components (e.g. the loop direction `{U0,U1,U2}` with `U0²+U1²+U2²=1`, or a
  ///        `{cos,sin}` pair). This is the general form of `sin²→1−cos²` and is the crux of the numeric
  ///        backend's compactness: it (a) collapses the bare-loop denominator `Σ(l1·Uμ)² = l1²·ΣUμ²` to
  ///        the **monomial `l1²`** so @ref divThroughMonomialAtoms cancels it (like inv's `rel`), and
  ///        (b) collapses the `U·U` factors transverse projectors generate. Value-preserving;
  ///        terminates because each rewrite strictly lowers the last component's exponent.
  inline MPoly reduce_units(const MPoly &p, const std::vector<std::vector<int>> &groups)
  {
    if (groups.empty()) return p;
    std::vector<std::pair<Mono, Cx>> out;
    std::vector<std::tuple<std::vector<int>, std::vector<int>, Cx>> work;
    for (const auto &[m, c] : p.t)
      work.push_back({m.e, m.atoms, c});
    while (!work.empty()) {
      auto [e, atoms, c] = std::move(work.back());
      work.pop_back();
      // find a unit group whose LAST component still has power >= 2 (the rewrite target)
      int groupIdx = -1;
      for (int gIdx = 0; gIdx < (int)groups.size(); ++gIdx)
        if (!groups[gIdx].empty() && e[groups[gIdx].back()] >= 2) {
          groupIdx = gIdx;
          break;
        }
      if (groupIdx < 0) {
        out.push_back({Mono{std::move(e), std::move(atoms)}, c});
        continue;
      }
      const std::vector<int> &group = groups[groupIdx];
      const int last = group.back();
      std::vector<int> base = e;
      base[last] -= 2;                                     // U_last^2 -> 1 - Σ_{μ<last} Uμ^2
      work.push_back({base, atoms, c});                    // the "+1" branch
      for (std::size_t i = 0; i + 1 < group.size(); ++i) { // the "-Uμ^2" branches
        std::vector<int> shifted = base;
        shifted[group[i]] += 2;
        work.push_back({std::move(shifted), atoms, Cx{-c.re, -c.im}});
      }
    }
    return MPoly::from_scratch(p.nsym, std::move(out));
  }

  /// @brief Numeric evaluation (validation only). `x[i]` = user symbol i; `atomVal[aid]` = value of
  ///        `1/D_aid` (the caller supplies the reciprocal already evaluated).
  inline Cx eval(const MPoly &p, const std::vector<double> &x, const std::vector<double> &atomVal)
  {
    Cx s{0, 0};
    for (const auto &[m, c] : p.t) {
      double mon = 1.0;
      for (int k = 0; k < p.nsym; ++k)
        for (int j = 0; j < m.e[k]; ++j)
          mon *= x[k];
      for (int aid : m.atoms)
        mon *= atomVal[aid];
      s = s + Cx{c.re * mon, c.im * mon};
    }
    return s;
  }

} // namespace numtracer::numeric
