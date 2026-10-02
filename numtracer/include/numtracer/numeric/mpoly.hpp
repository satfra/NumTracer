/// @file mpoly.hpp
/// @brief Multivariate polynomial over an arbitrary user-symbol set, carrying surviving inverse
///        (`1/k²`) atoms in the monomial key. The arithmetic core of the **numeric contraction
///        backend**: contract a diagram as 4×4 spinor matrix products whose entries are
///        these polynomials, so γ/metric/projector stay numeric and only the user's symbolic
///        momentum data enters the result — bounded by the matrix structure, not the `(2n−1)!!`
///        sp-monomial blowup of the symbolic `reduce_product` path.
///
/// A `MPoly` is a map from monomial → complex coefficient. A monomial is
///   (exponent vector over the `nsym` user symbols, a sorted multiset of inverse-atom ids).
/// The atom ids ride along so a transverse projector's `INV(k) = 1/k²` factor is tracked exactly, and
/// two passes cancel it against the numerator:
///   - @ref divThroughMonomialAtoms — denominator is a single **monomial** (e.g. a unit-direction loop
///     `k²=l1²`): cancels term-by-term wherever the numerator exponents dominate. The general
///     bare-loop cancellation, not a special case.
///   - @ref divThroughPolyAtoms — denominator is a genuine **polynomial** (a shifted line `k=l−q`,
///     `k²` non-monomial): terms are grouped by their atom multiset and the denominator is
///     trial-divided into the group, cancelling on a vanishing remainder. Dirac-trace numerators
///     routinely carry a factor of the very `k²` beneath them, so this fires often.
/// Only an atom that survives both is lowered to an `inv` env slot (a runtime division).
///
/// Map-based (sizes are tens of monomials with frame inputs — the regime where this backend wins);
/// the symbol *meaning* lives only in the kernel's `fill`, so the engine is frame-agnostic.
///
/// This header is the polynomial type alone; the 4×4 spinor matrix of `MPoly` (`Mat4`) and the
/// γ/slash builders that consume it live in `numeric/spinor_mat.hpp`.
#pragma once

#include "numtracer/core/config.hpp" // NT_THROW (exception-optional guard for -fno-exceptions builds)
#include "numtracer/core/cx.hpp"
#include "numtracer/numeric/stats.hpp"
#include "numtracer/third_party/gch/small_vector.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>
#include <map>
#include <cmath>

namespace numtracer::numeric
{

  using numtracer::Cx;

  /// Inline capacities for @ref Mono. `e` is ALWAYS exactly `nsym` long and `nsym` is the flow's env
  /// size (single digits to low tens), while `atoms` is usually empty and rarely more than a few — so
  /// with these buffers a monomial holds its whole payload inline and never touches the heap.
  ///
  /// This matters because `operator*` materialises all n·m product monomials: with `std::vector`
  /// members that was one heap allocation per product term, and profiling the generator run showed
  /// malloc/free at ~17% of total time (plus the pointer-chasing it forces on every sort comparison).
  /// Larger flows overflow these buffers gracefully — small_vector just falls back to the heap, i.e.
  /// exactly the old behaviour, so correctness never depends on the capacity being big enough.
  ///
  /// Exponents and atom ids are `int16_t`, not `int`. Both are small by construction — an exponent is
  /// a polynomial degree (single digits in practice) and an atom id indexes the flow's atom table — so
  /// 32 bits each is pure width. It matters because `std::sort` in @ref MPoly::from_scratch physically
  /// shuffles `pair<Mono, Cx>` objects and the monomial is the bulk of one: halving the exponent width
  /// halves what every swap moves and doubles how many monomials fit in cache.
  using MonoExpT = std::int16_t;
  using MonoAtomT = std::int16_t;

  inline constexpr unsigned kMonoAtomInline = 8;

  /// A stateless allocator identical to `std::allocator<T>` but advertising a 32-bit `size_type`.
  /// `gch::small_vector` stores its size and capacity in the allocator's `size_type` (packed to the
  /// smallest int that fits it — see small_vector.hpp), so this halves the @ref MonoAtoms header from
  /// 24 B (ptr + 8-byte size + 8-byte cap) to 16 B — one step of the `Mono` 112→96→56 B chain
  /// documented on @ref MonoExp. A monomial's atom count is single digits by construction, so the 2³²
  /// capacity ceiling is never approached; only the WIDTH of the bookkeeping changes, never a value —
  /// so every emitted kernel is bit-identical. Empty (EBO) like `std::allocator`, so it adds no bytes
  /// of its own.
  template <class T> struct NarrowAlloc : std::allocator<T> {
    using size_type = std::uint32_t;
    using difference_type = std::int32_t;
    template <class U> struct rebind { using other = NarrowAlloc<U>; };
    NarrowAlloc() = default;
    template <class U> NarrowAlloc(const NarrowAlloc<U> &) noexcept {}
  };

  /// @brief Packed exponent vector: the per-symbol exponents of a monomial stored as a fixed 128-bit
  ///        key (`kExpBits` = 5 bits/symbol → degree ≤ 31, up to 24 symbols) instead of the 72-byte
  ///        `small_vector`. This is the single largest RAM lever, and it makes `Mono::operator<`/`==`
  ///        two 64-bit integer compares instead of an nsym-long element walk.
  ///
  /// Size chain, measured: the stored polynomial term `pair<Mono,Cx>` went 128 B → 112 B (@ref
  /// NarrowAlloc on the atom list) → **72 B** here, i.e. −43.75% overall, which is why every
  /// `poly_bytes` estimate has a lower floor than it used to. The struct itself is 2 words + one
  /// pointer = 24 B, taking `Mono` to 56 B.
  ///
  /// Packing is **big-endian within each word** (symbol 0 in the high bits, 12 symbols per 64-bit word,
  /// so no symbol straddles the word boundary). That is deliberate: it makes the array's own ordering
  /// (`packed[0]` then `packed[1]`, integer compares) reproduce the OLD element-wise lexicographic order
  /// of the exponent vector EXACTLY — symbol 0 dominates, then symbol 1, … — for two monomials of equal
  /// nsym (always the case in one polynomial), with unused high slots zero in both. So sort/merge/lookup
  /// keep the identical order and the emitted kernel is byte-for-byte unchanged; only the STORAGE of the
  /// exponents changes, never a value or an ordering.
  ///
  /// The common inline range is nsym ≤ 24 with degree ≤ 31 — which covers every flow in practice (the
  /// same 24 the old `small_vector` inlined; the largest `nsym` across all committed flows is 6). It is
  /// NOT a hard cap: a symbol index ≥ 24, an exponent > 31, or a negative exponent transparently falls
  /// back to a heap `std::vector<MonoExpT>` holding the full exponent list, so nsym and degree are
  /// bounded only by `MonoExpT` (`int16_t`, ≤ 32767) — exactly like the old heap-expanding
  /// `small_vector`. The overflow path never triggers for the practical flows, so they keep the fast
  /// inline path; only an atypically large flow pays the indirection.
  struct MonoExp {
    static constexpr unsigned kExpBits = 5;
    static constexpr unsigned kPerWord = 64u / kExpBits;         ///< 12 symbols per 64-bit word (no straddle)
    static constexpr int kInlineSyms = 2 * static_cast<int>(kPerWord); ///< 24 symbols inline
    static constexpr MonoExpT kMaxExp = static_cast<MonoExpT>((1u << kExpBits) - 1u); ///< 31

    std::array<std::uint64_t, 2> packed{};             ///< packed inline exponents (big-endian per word)
    std::unique_ptr<std::vector<MonoExpT>> overflow;   ///< full exponent list; null ⇒ inline

    static unsigned wordOf(int k) { return static_cast<unsigned>(k) / kPerWord; }
    /// Big-endian slot: symbol 0 → shift 59 (high bits), symbol 11 → shift 4. Low nibble unused, so the
    /// integer compare of `packed[0]` then `packed[1]` is exactly the lexicographic order of symbols
    /// 0,1,2,….
    static unsigned shiftOf(int k) { return 64u - (static_cast<unsigned>(k) % kPerWord + 1u) * kExpBits; }
    MonoExpT unpackInline(int k) const
    {
      return static_cast<MonoExpT>((packed[wordOf(k)] >> shiftOf(k)) & static_cast<std::uint64_t>(kMaxExp));
    }
    /// Move all inline symbols into a fresh heap list (called when a set() first leaves the inline
    /// range). Afterwards `overflow` is authoritative and holds the 24 previously-inline exponents.
    void spillToHeap()
    {
      overflow = std::make_unique<std::vector<MonoExpT>>();
      overflow->reserve(static_cast<std::size_t>(kInlineSyms));
      for (int j = 0; j < kInlineSyms; ++j) overflow->push_back(unpackInline(j));
    }

    MonoExp() = default;
    /// `MonoExp(n, 0)`: the all-zero exponent (stays inline — an all-zero monomial never overflows,
    /// and a heap peer of length n>24 compares correctly against it since get() returns 0 past 24).
    MonoExp(int /*n*/, int fill) { assert(fill == 0); (void)fill; }
    /// Pack from an iterator range of integer exponents (the generated component table hands a vector).
    template <class It> MonoExp(It b, It e)
    {
      int k = 0;
      for (; b != e; ++b, ++k) set(k, static_cast<MonoExpT>(*b));
    }

    // `overflow` owns its list, so the moves are the compiler's; only the deep COPY is hand-written.
    MonoExp(MonoExp &&) noexcept = default;
    MonoExp &operator=(MonoExp &&) noexcept = default;
    MonoExp(const MonoExp &o)
        : packed(o.packed),
          overflow(o.overflow ? std::make_unique<std::vector<MonoExpT>>(*o.overflow) : nullptr)
    {
    }
    MonoExp &operator=(const MonoExp &o)
    {
      if (this != &o) {
        auto copy = o.overflow ? std::make_unique<std::vector<MonoExpT>>(*o.overflow) : nullptr;
        packed = o.packed;
        overflow = std::move(copy);
      }
      return *this;
    }

    /// Number of symbol slots the slow (heap-involved) compares must walk: the heap list's length, or
    /// the inline ceiling. Symbols past a monomial's own extent read 0, so an over-long bound is
    /// harmless.
    int symbolCount() const { return overflow ? static_cast<int>(overflow->size()) : kInlineSyms; }

    MonoExpT get(int k) const
    {
      if (overflow)
        return k < static_cast<int>(overflow->size()) ? (*overflow)[static_cast<std::size_t>(k)] : MonoExpT(0);
      return k < kInlineSyms ? unpackInline(k) : MonoExpT(0);
    }
    void set(int k, MonoExpT v)
    {
      assert(k >= 0); // the intended domain; a negative exponent is still stored faithfully below
      if (overflow) {
        if (k >= static_cast<int>(overflow->size())) overflow->resize(static_cast<std::size_t>(k) + 1, 0);
        (*overflow)[static_cast<std::size_t>(k)] = v;
        return;
      }
      if (k < kInlineSyms && v >= 0 && v <= kMaxExp) { // fast inline path
        const unsigned sh = shiftOf(k);
        const std::uint64_t mask = static_cast<std::uint64_t>(kMaxExp) << sh;
        packed[wordOf(k)] = (packed[wordOf(k)] & ~mask) | ((static_cast<std::uint64_t>(v) & kMaxExp) << sh);
        return;
      }
      // A zero past the inline range is already what get() reports, so storing it changes nothing —
      // and returning here is what keeps `operator*` allocation-free on an nsym > kInlineSyms flow,
      // where its `for (k < ns) m.e[k] = …` writes a zero to every high slot of every product monomial.
      if (v == 0) return;
      // Genuine overflow: symbol index ≥ kInlineSyms, exponent > kMaxExp, or NEGATIVE. The last case
      // matters — masking a negative into 5 bits would silently store a bogus positive exponent, so it
      // takes the faithful heap representation instead. (Both subtracting call sites,
      // divThroughMonomialAtoms and divThroughPolyAtoms, guard against going negative; this is the
      // backstop for when that invariant moves.)
      spillToHeap();
      if (k >= static_cast<int>(overflow->size())) overflow->resize(static_cast<std::size_t>(k) + 1, 0);
      (*overflow)[static_cast<std::size_t>(k)] = v;
    }

    /// Mutable element proxy, so existing `e[k] = / += / -=` call sites keep working verbatim.
    struct Ref {
      MonoExp *owner;
      int slot;
      operator MonoExpT() const { return owner->get(slot); }
      Ref &operator=(MonoExpT v) { owner->set(slot, v); return *this; }
      Ref &operator+=(MonoExpT v) { owner->set(slot, static_cast<MonoExpT>(owner->get(slot) + v)); return *this; }
      Ref &operator-=(MonoExpT v) { owner->set(slot, static_cast<MonoExpT>(owner->get(slot) - v)); return *this; }
    };
    MonoExpT operator[](int k) const { return get(k); }
    Ref operator[](int k) { return Ref{this, k}; }

    /// Mask of every field's TOP bit (bit `shiftOf(k)+kExpBits-1`). A 64-bit add of two packed
    /// words is the exact fieldwise sum IFF no field's sum carries out of its top bit; the classic
    /// SWAR carry-out detect `(a&b) | ((a|b) & ~sum)` tested against this mask decides that. Both
    /// inputs are ≤ 31 per field, so a field sum is ≤ 62 and at most one carry bit exists.
    static constexpr std::uint64_t kFieldTop = [] {
      std::uint64_t m = 0;
      for (unsigned k = 0; k < kPerWord; ++k)
        m |= 1ULL << (64u - (k + 1u) * kExpBits + (kExpBits - 1u));
      return m;
    }();

    bool operator<(const MonoExp &o) const
    {
      // fast path: both inline (symbols ≥ kInlineSyms are 0 in both)
      if (!overflow && !o.overflow) return packed < o.packed;
      const int n = std::max(symbolCount(), o.symbolCount());
      for (int k = 0; k < n; ++k) {
        const MonoExpT a = get(k), b = o.get(k);
        if (a != b) return a < b;
      }
      return false;
    }
    bool operator==(const MonoExp &o) const
    {
      if (!overflow && !o.overflow) return packed == o.packed;
      const int n = std::max(symbolCount(), o.symbolCount());
      for (int k = 0; k < n; ++k)
        if (get(k) != o.get(k)) return false;
      return true;
    }
    bool operator!=(const MonoExp &o) const { return !(*this == o); }
  };

  using MonoAtoms = gch::small_vector<MonoAtomT, kMonoAtomInline, NarrowAlloc<MonoAtomT>>;

  struct Mono;
  /// Scratch list of (monomial, coeff) handed to @ref MPoly::from_scratch. Measured on a real flow,
  /// from_scratch is called tens of millions of times with a MEAN of ~4 terms — it is not one big
  /// sort, it is a flood of tiny ones — so a heap-backed scratch vector paid an allocation per call.
  /// Inline storage covers the common case; larger products fall back to the heap as before.
  inline constexpr unsigned kMPolyScratchInline = 8;

  /// @brief A monomial: exponents over the user symbols plus a sorted multiset of inverse-atom ids.
  struct Mono {
    MonoExp e;       ///< length nsym, exponent of each user symbol
    MonoAtoms atoms; ///< sorted (with multiplicity) ids of surviving `1/k²` atoms
    bool operator<(const Mono &o) const
    {
      if (e != o.e) return e < o.e;
      return atoms < o.atoms;
    }
    bool operator==(const Mono &o) const { return e == o.e && atoms == o.atoms; }
  };

  using MPolyScratch = gch::small_vector<std::pair<Mono, Cx>, kMPolyScratchInline>;

  // The `nsym`-carrying construction API is closed behind these friends (see @ref LorentzEnv): every
  // polynomial in one trace must share an `nsym`, so the factories that bake it in are private and the
  // ONLY sanctioned way to mint one is a @ref LorentzEnv (which holds a fixed `nsym`). @ref MPolyFactory
  // is a tiny internal attorney that re-exposes the same factories to the trusted cross-header engine
  // code (contraction / trace-fold), which already threads a single `nsym` and must not route through a
  // user-facing env. Only names are needed here; the definitions live below / in `numeric/env.hpp`.
  class LorentzEnv;
  struct MPolyFactory;

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
    /// Heap-backed on purpose. Giving `t` inline storage (capacity 4) was measured: it bought ~2% of
    /// run time but pushed sizeof(Mat4) from ~0.5 KB to ~8.7 KB and peak RSS up ~20% — a bad trade on
    /// the dense flows, which are memory-bound before they are time-bound.
    std::vector<std::pair<Mono, Cx>> terms; ///< sorted by Mono, like terms combined, no zeros

    // Sanctioned construction paths (see the note above @ref MPoly). The env and the internal attorney
    // reach the nsym-taking ctor/factories; the in-header arithmetic operators keep constructing result
    // polynomials directly (they already have a definite `nsym` from their operands), so their bodies
    // stay byte-for-byte unchanged — the hot path is untouched.
    friend class LorentzEnv;
    friend struct MPolyFactory;
    friend MPoly operator+(const MPoly &a, const MPoly &b);
    friend MPoly operator-(const MPoly &a, const MPoly &b);
    friend MPoly operator*(const MPoly &a, const MPoly &b);
    friend MPoly divThroughMonomialAtoms(const MPoly &p, const std::vector<MPoly> &atomDen);
    friend MPoly reduce_units(const MPoly &p, const std::vector<std::vector<int>> &groups);
    friend MPoly divThroughPolyAtoms(const MPoly &p, const std::vector<MPoly> &atomDen);
    // The rebuild halves of the two reductions above, split out so their `const&` and `&&` overloads
    // share one body. Same trust level as the functions they were extracted from — the `&&` overloads
    // themselves need no friendship, since they only move and read public members.
    friend MPoly dmaRebuild(const MPoly &p, const std::vector<MPoly> &atomDen);
    friend MPoly reduceUnitsRebuild(const MPoly &p, const std::vector<std::vector<int>> &groups);

    MPoly() = default;

  private:
    // Bare-`nsym` construction — reachable only through @ref LorentzEnv / @ref MPolyFactory (friends).
    // Making these private is what turns "every MPoly in one trace shares an nsym" from a debug-time
    // assert into a compile-time guarantee: outside the sanctioned env you cannot mint a non-empty MPoly
    // with a hand-picked nsym. The empty default ctor above stays public (an nsym==0 zero used by the
    // operator short-circuits and by std::array/std::vector default members).
    explicit MPoly(int ns) : nsym(ns) {}

    /// Build from an unsorted scratch list of (monomial, coeff): sort then combine adjacent equals.
    ///
    /// The sort runs `std::sort` over 24-byte proxies (packed exponents + scratch index), not over the
    /// 72-byte terms themselves. `std::sort`'s permutation depends only on its comparison results, and
    /// the proxy comparator reproduces `Mono::operator<` exactly, so the permutation — and with it the
    /// order like terms are summed in — is identical to sorting `s` directly. A scratch that is
    /// already strictly increasing (e.g. a product by a constant) skips the sort: with distinct keys
    /// every sort is the identity.
    static MPoly from_scratch(int ns, MPolyScratch s)
    {
      NT_STAT_ADD(fs_calls, 1);
      NT_STAT_ADD(fs_terms_in, s.size());
      MPoly p(ns);
      p.terms.reserve(s.size());
      auto absorb = [&p](std::pair<Mono, Cx> &kv) {
        if (kv.second.re == 0 && kv.second.im == 0) return;
        if (!p.terms.empty() && p.terms.back().first == kv.first) {
          p.terms.back().second = p.terms.back().second + kv.second;
          if (p.terms.back().second.re == 0 && p.terms.back().second.im == 0) p.terms.pop_back();
        } else {
          p.terms.push_back(std::move(kv));
        }
      };
      bool sorted = true;
      for (std::size_t i = 1; sorted && i < s.size(); ++i)
        sorted = s[i - 1].first < s[i].first;
      if (sorted) {
        for (auto &kv : s) absorb(kv);
        return p;
      }
      struct Key {
        std::uint64_t p0, p1;
        std::uint32_t idx;
        bool inlineExp; ///< false ⇒ heap exponents, compare through the full Mono
      };
      gch::small_vector<Key, 16> keys;
      keys.reserve(s.size());
      for (std::size_t i = 0; i < s.size(); ++i) {
        const MonoExp &e = s[i].first.e;
        keys.push_back({e.packed[0], e.packed[1], static_cast<std::uint32_t>(i), !e.overflow});
      }
      std::sort(keys.begin(), keys.end(), [&s](const Key &a, const Key &b) {
        const Mono &ma = s[a.idx].first, &mb = s[b.idx].first;
        if (!(a.inlineExp && b.inlineExp)) return ma < mb;
        if (a.p0 != b.p0) return a.p0 < b.p0;
        if (a.p1 != b.p1) return a.p1 < b.p1;
        return ma.atoms < mb.atoms;
      });
      for (const Key &k : keys) absorb(s[k.idx]);
      return p;
    }

    /// True when no term carries inverse atoms or heap-stored exponents — the case @ref operator*
    /// handles with packed keys alone.
    static bool atomFreeInline(const MPoly &p)
    {
      for (const auto &kv : p.terms)
        if (kv.first.e.overflow || !kv.first.atoms.empty()) return false;
      return true;
    }

    /// One term of an atom-free product in key form: packed exponents and coefficient.
    struct KeyedTerm {
      std::uint64_t p0, p1;
      Cx c;
    };
    using KeyedTerms = gch::small_vector<KeyedTerm, 16>;

    /// The product `a·b` of atom-free, inline-exponent operands as sorted, like-terms-combined
    /// @ref KeyedTerm s — exactly the terms @ref operator* returns, without building any Mono.
    /// Returns false on a field carry (an exponent > 31), which needs the heap representation.
    static bool keyedAtomFreeTerms(const MPoly &a, const MPoly &b, KeyedTerms &out)
    {
      const std::size_t na = a.terms.size(), nb = b.terms.size();
      struct Key {
        std::uint64_t p0, p1;
        std::uint32_t idx;
      };
      gch::small_vector<Key, 16> keys;
      gch::small_vector<Cx, 16> coeff;
      keys.reserve(na * nb);
      coeff.reserve(na * nb);
      bool carry = false;
      for (const auto &[ma, ca] : a.terms) {
        for (const auto &[mb, cb] : b.terms) {
          std::uint64_t sum[2];
          for (int word = 0; word < 2; ++word) {
            const std::uint64_t wa = ma.e.packed[word], wb = mb.e.packed[word];
            sum[word] = wa + wb;
            carry |= (((wa & wb) | ((wa | wb) & ~sum[word])) & MonoExp::kFieldTop) != 0;
          }
          keys.push_back({sum[0], sum[1], static_cast<std::uint32_t>(keys.size())});
          coeff.push_back(ca * cb);
        }
        if (carry) return false;
      }
      NT_STAT_ADD(mul_keyed, 1);
      auto less = [](const Key &x, const Key &y) { return x.p0 != y.p0 ? x.p0 < y.p0 : x.p1 < y.p1; };
      bool sorted = true;
      for (std::size_t k = 1; sorted && k < keys.size(); ++k)
        sorted = less(keys[k - 1], keys[k]);
      if (!sorted) std::sort(keys.begin(), keys.end(), less);
      out.clear();
      out.reserve(keys.size());
      for (const Key &k : keys) {
        const Cx c = coeff[k.idx];
        if (c.re == 0 && c.im == 0) continue;
        if (!out.empty() && out.back().p0 == k.p0 && out.back().p1 == k.p1) {
          out.back().c = out.back().c + c;
          if (out.back().c.re == 0 && out.back().c.im == 0) out.pop_back();
        } else {
          out.push_back({k.p0, k.p1, c});
        }
      }
      return true;
    }

    /// Build the polynomial from key-form terms (atom-free, inline exponents).
    static MPoly fromKeyed(int ns, const KeyedTerm *t, std::size_t n)
    {
      MPoly p(ns);
      p.terms.resize(n);
      for (std::size_t k = 0; k < n; ++k) {
        p.terms[k].first.e.packed = {t[k].p0, t[k].p1};
        p.terms[k].second = t[k].c;
      }
      return p;
    }

    /// Keyed product for atom-free operands (see @ref operator*); `nullopt` on a field carry.
    [[gnu::noinline]] static std::optional<MPoly> mulKeyedAtomFree(const MPoly &a, const MPoly &b)
    {
      const std::size_t na = a.terms.size(), nb = b.terms.size();
      struct Key {
        std::uint64_t p0, p1;
        std::uint32_t idx;
      };
      gch::small_vector<Key, 16> keys;
      gch::small_vector<Cx, 16> coeff;
      keys.reserve(na * nb);
      coeff.reserve(na * nb);
      bool carry = false;
      for (const auto &[ma, ca] : a.terms) {
        for (const auto &[mb, cb] : b.terms) {
          std::uint64_t sum[2];
          for (int word = 0; word < 2; ++word) {
            const std::uint64_t wa = ma.e.packed[word], wb = mb.e.packed[word];
            sum[word] = wa + wb;
            carry |= (((wa & wb) | ((wa | wb) & ~sum[word])) & MonoExp::kFieldTop) != 0;
          }
          keys.push_back({sum[0], sum[1], static_cast<std::uint32_t>(keys.size())});
          coeff.push_back(ca * cb);
        }
        if (carry) return std::nullopt;
      }
      NT_STAT_ADD(mul_keyed, 1);
      auto less = [](const Key &x, const Key &y) { return x.p0 != y.p0 ? x.p0 < y.p0 : x.p1 < y.p1; };
      bool sorted = true;
      for (std::size_t k = 1; sorted && k < keys.size(); ++k)
        sorted = less(keys[k - 1], keys[k]);
      if (!sorted) std::sort(keys.begin(), keys.end(), less);
      MPoly p(a.nsym);
      p.terms.reserve(keys.size());
      std::uint64_t back0 = 0, back1 = 0;
      for (const Key &k : keys) {
        const Cx c = coeff[k.idx];
        if (c.re == 0 && c.im == 0) continue;
        if (!p.terms.empty() && back0 == k.p0 && back1 == k.p1) {
          Cx &acc = p.terms.back().second;
          acc = acc + c;
          if (acc.re == 0 && acc.im == 0) {
            p.terms.pop_back();
            if (!p.terms.empty()) {
              back0 = p.terms.back().first.e.packed[0];
              back1 = p.terms.back().first.e.packed[1];
            }
          }
        } else {
          p.terms.emplace_back();
          p.terms.back().first.e.packed = {k.p0, k.p1};
          p.terms.back().second = c;
          back0 = k.p0;
          back1 = k.p1;
        }
      }
      return p;
    }

    /// `x0·y0 + x1·y1`, bit-identical to evaluating it with @ref operator* and @ref operator+, for
    /// the Dirac fold's 2×2 block products. With atom-free operands the two products stay in key form
    /// and are merged there, so only the final terms are built as Monos.
    [[gnu::noinline]] static MPoly mulAdd(const MPoly &x0, const MPoly &y0, const MPoly &x1, const MPoly &y1);

    /// Keyed product for operands that carry atoms: the distinct merged atom lists are built once and
    /// ranked in MonoAtoms order, and the key is (packed exponent sum, atom rank), which compares exactly
    /// like the full monomial. `nullopt` ⇒ use the scratch path (heap exponents, a field carry, or more
    /// than kMaxLists distinct atom lists in one operand).
    [[gnu::noinline]] static std::optional<MPoly> mulKeyed(const MPoly &a, const MPoly &b)
    {
      constexpr std::size_t kMaxLists = 64;
      // Distinct atom lists of one operand, and each term's index into them.
      auto atomLists = [](const MPoly &p, gch::small_vector<const MonoAtoms *, 8> &lists,
                          gch::small_vector<std::uint32_t, 16> &of) {
        of.reserve(p.terms.size());
        for (const auto &kv : p.terms) {
          if (kv.first.e.overflow) return false;
          std::size_t x = 0;
          while (x < lists.size() && !(*lists[x] == kv.first.atoms)) ++x;
          if (x == lists.size()) {
            if (lists.size() == kMaxLists) return false;
            lists.push_back(&kv.first.atoms);
          }
          of.push_back(static_cast<std::uint32_t>(x));
        }
        return true;
      };
      gch::small_vector<const MonoAtoms *, 8> la, lb;
      gch::small_vector<std::uint32_t, 16> ofA, ofB;
      if (!atomLists(a, la, ofA) || !atomLists(b, lb, ofB)) return std::nullopt;
      // Merged multiset for every list pair, ranked by MonoAtoms order (equal lists share a rank).
      const std::size_t nPairs = la.size() * lb.size();
      gch::small_vector<MonoAtoms, 1> merged(nPairs);
      for (std::size_t x = 0; x < la.size(); ++x)
        for (std::size_t y = 0; y < lb.size(); ++y) {
          const MonoAtoms &u = *la[x], &v = *lb[y];
          MonoAtoms &m = merged[x * lb.size() + y];
          std::size_t i = 0, j = 0;
          while (i < u.size() && j < v.size()) m.push_back(u[i] <= v[j] ? u[i++] : v[j++]);
          while (i < u.size()) m.push_back(u[i++]);
          while (j < v.size()) m.push_back(v[j++]);
        }
      gch::small_vector<std::uint32_t, 1> rank(nPairs, 0);
      if (nPairs > 1) {
        gch::small_vector<std::uint32_t, 16> order(nPairs);
        for (std::size_t q = 0; q < nPairs; ++q) order[q] = static_cast<std::uint32_t>(q);
        std::sort(order.begin(), order.end(), [&](std::uint32_t x, std::uint32_t y) { return merged[x] < merged[y]; });
        std::uint32_t r = 0;
        for (std::size_t q = 0; q < nPairs; ++q) {
          if (q > 0 && merged[order[q - 1]] < merged[order[q]]) ++r;
          rank[order[q]] = r;
        }
      }
      gch::small_vector<std::uint32_t, 1> repOfRank(nPairs, 0); // one merged list per rank
      for (std::size_t q = 0; q < nPairs; ++q) repOfRank[rank[q]] = static_cast<std::uint32_t>(q);

      struct Key {
        std::uint64_t p0, p1;
        std::uint32_t rank, idx;
      };
      const std::size_t na = a.terms.size(), nb = b.terms.size();
      gch::small_vector<Key, 16> keys;
      gch::small_vector<Cx, 16> coeff;
      keys.reserve(na * nb);
      coeff.reserve(na * nb);
      for (std::size_t i = 0; i < na; ++i) {
        const auto &[ma, ca] = a.terms[i];
        const std::size_t rowA = ofA[i] * lb.size();
        bool carry = false;
        for (std::size_t j = 0; j < nb; ++j) {
          const auto &[mb, cb] = b.terms[j];
          std::uint64_t sum[2];
          for (int word = 0; word < 2; ++word) {
            const std::uint64_t wa = ma.e.packed[word], wb = mb.e.packed[word];
            sum[word] = wa + wb;
            carry |= (((wa & wb) | ((wa | wb) & ~sum[word])) & MonoExp::kFieldTop) != 0;
          }
          keys.push_back({sum[0], sum[1], rank[rowA + ofB[j]], static_cast<std::uint32_t>(keys.size())});
          coeff.push_back(ca * cb);
        }
        if (carry) return std::nullopt;
      }
      NT_STAT_ADD(mul_keyed, 1);
      auto less = [](const Key &x, const Key &y) {
        if (x.p0 != y.p0) return x.p0 < y.p0;
        if (x.p1 != y.p1) return x.p1 < y.p1;
        return x.rank < y.rank;
      };
      bool sorted = true;
      for (std::size_t k = 1; sorted && k < keys.size(); ++k)
        sorted = less(keys[k - 1], keys[k]);
      if (!sorted) std::sort(keys.begin(), keys.end(), less);
      MPoly p(a.nsym);
      p.terms.reserve(keys.size());
      gch::small_vector<const Key *, 16> backKey; // the key of each result term, for the like-term test
      for (const Key &k : keys) {
        const Cx c = coeff[k.idx];
        if (c.re == 0 && c.im == 0) continue;
        if (!backKey.empty() && backKey.back()->p0 == k.p0 && backKey.back()->p1 == k.p1 &&
            backKey.back()->rank == k.rank) {
          Cx &acc = p.terms.back().second;
          acc = acc + c;
          if (acc.re == 0 && acc.im == 0) {
            p.terms.pop_back();
            backKey.pop_back();
          }
          continue;
        }
        p.terms.emplace_back();
        p.terms.back().first.e.packed = {k.p0, k.p1};
        p.terms.back().first.atoms = merged[repOfRank[k.rank]];
        p.terms.back().second = c;
        backKey.push_back(&k);
      }
      return p;
    }

    static MPoly constant(int ns, Cx c)
    {
      MPoly p(ns);
      if (!(c.re == 0 && c.im == 0)) p.terms.push_back({Mono{MonoExp(ns, 0), {}}, c});
      return p;
    }

    /// Scale every coefficient of `p` by the constant `c`. This is exactly `constant(ns, c) * p`, but
    /// without the n·m scratch and the `std::sort` that `operator*` pays: multiplying by a constant
    /// leaves every monomial (hence `p`'s existing sort order and its distinctness) unchanged, so there
    /// is nothing to sort and no like terms to collect. `c * coeff` matches the operand order of the
    /// `constant(ns,c) * p` it replaces — and `Cx` multiply is componentwise, so `c*coeff == coeff*c`
    /// bit-for-bit, making this bit-identical to BOTH the `constant * p` and `p * constant` call sites.
    /// A nonzero `c` times a stored (nonzero) coeff is nonzero, so no term can vanish; `c == 0` yields
    /// the empty (zero) polynomial, exactly as `constant(ns, 0) * p` does.
    static MPoly scaled(int ns, const MPoly &p, Cx c)
    {
      MPoly r(ns);
      if (c.re == 0 && c.im == 0) return r;
      r.terms.reserve(p.terms.size());
      for (const auto &kv : p.terms)
        r.terms.push_back({kv.first, c * kv.second});
      return r;
    }
    /// The i-th user symbol (coefficient 1).
    static MPoly var(int ns, int i)
    {
      MPoly p(ns);
      MonoExp e(ns, 0);
      e[i] = 1;
      p.terms.push_back({Mono{std::move(e), {}}, Cx{1, 0}});
      return p;
    }
    /// A single monomial `c · ∏ x_k^{e_k}` (no inverse atoms). Used by the generated component table.
    /// Takes a `std::vector` because the GENERATED component table hands one over as a braced list;
    /// it is converted into the monomial's inline storage here (once per table entry, not on any hot
    /// path).
    static MPoly mono(int ns, const std::vector<int> &e, Cx c)
    {
      MPoly p(ns);
      if (!(c.re == 0 && c.im == 0)) p.terms.push_back({Mono{MonoExp(e.begin(), e.end()), {}}, c});
      return p;
    }
    /// A bare inverse atom `1/D` (atom id `aid`), coefficient 1.
    ///
    /// The narrowing to @ref MonoAtomT (`int16_t`) is the ONE place an atom id enters the monomial
    /// key, and it was unguarded: `aid >= 32768` wraps to a different — possibly negative — id, so
    /// the term then carries somebody else's denominator and @ref divThroughMonomialAtoms cancels
    /// against the wrong `atomDen` entry. Silent, and value-wrong. The id space is small by
    /// construction (one per projector denominator in the flow), so tripping this means the front
    /// end changed, not that the bound is too tight.
    static MPoly atom(int ns, int aid)
    {
      if (aid < 0 || aid > static_cast<int>(std::numeric_limits<MonoAtomT>::max()))
        NT_THROW(std::runtime_error, "MPoly::atom: atom id out of MonoAtomT (int16) range — it would "
                                     "wrap silently and alias another denominator");
      MPoly p(ns);
      p.terms.push_back({Mono{MonoExp(ns, 0), MonoAtoms{static_cast<MonoAtomT>(aid)}}, Cx{1, 0}});
      return p;
    }

  public:
    int size() const { return terms.size(); }
    bool empty() const { return terms.empty(); }

    /// Insert/accumulate one term (keeps `terms` sorted). O(n) shift — used only on the incremental paths
    /// (component builders); the hot `operator*`/`operator+`/reductions build whole vectors at once.
    void addTerm(const Mono &m, Cx c)
    {
      if (c.re == 0 && c.im == 0) return;
      auto it = std::lower_bound(terms.begin(), terms.end(), m,
                                 [](const std::pair<Mono, Cx> &a, const Mono &k) { return a.first < k; });
      if (it != terms.end() && it->first == m) {
        it->second = it->second + c;
        if (it->second.re == 0 && it->second.im == 0) terms.erase(it);
      } else {
        terms.insert(it, {m, c});
      }
    }
  };

  /// @brief Internal attorney re-exposing the private @ref MPoly factories to the trusted cross-header
  ///        engine code (spinor matrices, contraction, trace-fold). Those functions already carry one
  ///        definite `nsym` and must not depend on the user-facing @ref LorentzEnv, but they cannot be
  ///        friended by name (templates / DiracNet-heavy signatures across headers), so the friend
  ///        surface is localised to this one struct. NOT part of the public API — call sites outside the
  ///        engine construct polynomials through @ref LorentzEnv.
  struct MPolyFactory {
    static MPoly zero(int ns) { return MPoly(ns); }
    static MPoly constant(int ns, Cx c) { return MPoly::constant(ns, c); }
    static MPoly scaled(int ns, const MPoly &p, Cx c) { return MPoly::scaled(ns, p, c); }
    static MPoly atom(int ns, int aid) { return MPoly::atom(ns, aid); }
    static MPoly from_scratch(int ns, MPolyScratch s) { return MPoly::from_scratch(ns, std::move(s)); }
    static MPoly mul_add(const MPoly &x0, const MPoly &y0, const MPoly &x1, const MPoly &y1)
    {
      return MPoly::mulAdd(x0, y0, x1, y1);
    }
  };

  inline MPoly operator+(const MPoly &a, const MPoly &b)
  {
    if (a.terms.empty()) {
      NT_STAT_ADD(add_empty, 1); // this return COPY-constructs the surviving polynomial
      return b;
    }
    if (b.terms.empty()) {
      NT_STAT_ADD(add_empty, 1);
      return a;
    }
    NT_STAT_ADD(add_calls, 1);
    // Both operands carry terms, so their symbol spaces must agree: the merge below walks the two
    // exponent vectors slot-for-slot, and a mismatch would read out of bounds (a default-constructed
    // zero with nsym == 0 is allowed — it is caught by the empty checks above). Debug-only; compiles
    // out under NDEBUG.
    assert(a.nsym == b.nsym);
    const int ns = a.nsym ? a.nsym : b.nsym;
    MPoly r(ns);
    r.terms.reserve(a.terms.size() + b.terms.size());
    std::size_t i = 0, j = 0;
    while (i < a.terms.size() && j < b.terms.size()) {
      if (a.terms[i].first < b.terms[j].first)
        r.terms.push_back(a.terms[i++]);
      else if (b.terms[j].first < a.terms[i].first)
        r.terms.push_back(b.terms[j++]);
      else {
        Cx s = a.terms[i].second + b.terms[j].second;
        if (!(s.re == 0 && s.im == 0)) r.terms.push_back({a.terms[i].first, s});
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
  /// Rvalue overloads of operator+: the hot paths add PRVALUE products (`x[0]*y[0] + x[1]*y[2]` in
  /// the Weyl-block fold) where one side is very often the zero polynomial — a γ Weyl block is
  /// diagonal-or-antidiagonal, so half its entry products are empty. The const& overload must COPY
  /// the surviving side in that case (measured: 70-85% of ALL operator+ calls); these move it.
  /// Value- and byte-identical: the merge path and term order are untouched.
  inline MPoly operator+(MPoly &&a, MPoly &&b)
  {
    if (a.terms.empty()) {
      NT_STAT_ADD(add_moved, 1);
      return std::move(b);
    }
    if (b.terms.empty()) {
      NT_STAT_ADD(add_moved, 1);
      return std::move(a);
    }
    return static_cast<const MPoly &>(a) + static_cast<const MPoly &>(b);
  }
  inline MPoly operator+(MPoly &&a, const MPoly &b)
  {
    if (b.terms.empty()) {
      NT_STAT_ADD(add_moved, 1);
      return std::move(a);
    }
    return static_cast<const MPoly &>(a) + b;
  }
  inline MPoly operator+(const MPoly &a, MPoly &&b)
  {
    if (a.terms.empty()) {
      NT_STAT_ADD(add_moved, 1);
      return std::move(b);
    }
    return a + static_cast<const MPoly &>(b);
  }
  inline MPoly operator-(const MPoly &a, const MPoly &b)
  {
    assert(a.terms.empty() || b.terms.empty() || a.nsym == b.nsym); // see operator+; debug-only
    MPoly nb(b.nsym);
    nb.terms.reserve(b.terms.size());
    for (const auto &kv : b.terms)
      nb.terms.push_back({kv.first, Cx{-kv.second.re, -kv.second.im}});
    return a + std::move(nb); // rvalue overload: an empty `a` MOVES the negation instead of copying it
  }
  /// Above this many product terms (`|a|·|b|`), `operator*` materialises its scratch in CHUNKS of the
  /// outer operand and folds the chunks, so peak scratch is bounded by ~`kMulMaxScratch` instead of the
  /// full `|a|·|b|`. This is the single-contraction RAM lever: `operator*` builds every `n·m` product
  /// monomial before `from_scratch` collapses like terms, so on a "collapsing" multiply (result ≪
  /// scratch) the transient peak — not the result — is what spikes RAM. Below the threshold the exact
  /// byte-for-byte path runs, so the mean-few-term hot path and every small/test flow are untouched.
  /// Blocking only REASSOCIATES the like-term sums across chunk boundaries (≤ 1 ulp, exactly as the
  /// already-shipped phase-B tree fold does); the monomial SET is identical.
  inline constexpr std::size_t kMulMaxScratch = std::size_t(1) << 20;
  /// Smallest product size (`|a|·|b|`) that takes the atom-ranked keyed multiply.
  inline constexpr std::size_t kMulKeyedMin = 32;

  inline MPoly operator*(const MPoly &a, const MPoly &b)
  {
    const int ns = a.nsym ? a.nsym : b.nsym;
    if (a.terms.empty() || b.terms.empty()) {
      NT_STAT_ADD(mul_empty, 1);
      return MPoly(ns);
    }
    assert(a.nsym == b.nsym); // both carry terms ⇒ symbol spaces must match; see operator+ (debug-only)
    const std::size_t na = a.terms.size(), nb = b.terms.size();
    NT_STAT_ADD(mul_calls, 1);
    NT_STAT_ADD(mul_prod_terms, na * nb);

    // Emit the product monomial ma·mb (coefficient ca·cb) into scratch `s`. Build IN PLACE (the old
    // code kept a scratch `e` and copied it in, an allocation per product term).
    auto emit = [ns](MPolyScratch &s, const Mono &ma, Cx ca, const Mono &mb, Cx cb) {
      auto &slot = s.emplace_back(Mono{}, ca * cb);
      Mono &m = slot.first;
      // Fieldwise exponent add. Fast path: both operands inline ⇒ `packed[w] + packed[w]` IS the
      // per-symbol sum whenever no 5-bit field carries out (SWAR carry-out test against the field
      // top-bit mask) — two 64-bit adds replace the nsym unpack/shift/mask/repack round trips of
      // the innermost loop of the innermost operation. Exact integer arithmetic: any carry (a
      // field sum > 31) falls back to the per-symbol path, which spills to the heap exactly as
      // before, so values and ordering are bit-identical in every case.
      bool fast = !ma.e.overflow && !mb.e.overflow;
      if (fast) {
        // `wa`/`wb`, not `a`/`b`: those name this operator's two POLYNOMIAL operands in the enclosing
        // scope, and these are one 64-bit word of the packed exponent vectors of one monomial each.
        for (int word = 0; word < 2; ++word) {
          const std::uint64_t wa = ma.e.packed[word], wb = mb.e.packed[word], sum = wa + wb;
          if (((wa & wb) | ((wa | wb) & ~sum)) & MonoExp::kFieldTop) {
            fast = false;
            m.e.packed = {0, 0};
            break;
          }
          m.e.packed[word] = sum;
        }
      }
      if (!fast)
        for (int k = 0; k < ns; ++k) m.e[k] = ma.e[k] + mb.e[k]; // packed exponents: set each directly
      // Only ask for capacity when the merge would overflow the inline buffer — most monomials carry
      // NO atoms, and an unconditional reserve() here cost ~7% of the run.
      const std::size_t nat = ma.atoms.size() + mb.atoms.size();
      if (nat > kMonoAtomInline) m.atoms.reserve(nat);
      std::size_t i = 0, j = 0; // merge the two sorted atom multisets
      while (i < ma.atoms.size() && j < mb.atoms.size())
        m.atoms.push_back(ma.atoms[i] <= mb.atoms[j] ? ma.atoms[i++] : mb.atoms[j++]);
      while (i < ma.atoms.size()) m.atoms.push_back(ma.atoms[i++]);
      while (j < mb.atoms.size()) m.atoms.push_back(mb.atoms[j++]);
    };

    // A factor that is a single constant term (a γ or C matrix entry in the Dirac fold) leaves every
    // monomial of the other factor unchanged, so the product is already sorted and distinct: copy the
    // terms, multiply the coefficients in the general path's operand order (ca·cb), and drop exact-zero
    // products as from_scratch would. Bit-identical to the scratch path, without the scratch.
    auto isConstant = [](const MPoly &p) {
      if (p.terms.size() != 1) return false;
      const Mono &m = p.terms[0].first;
      return !m.e.overflow && m.e.packed[0] == 0 && m.e.packed[1] == 0 && m.atoms.empty();
    };
    if (nb == 1 && isConstant(b)) {
      NT_STAT_ADD(mul_const, 1);
      MPoly r(ns);
      r.terms.reserve(na);
      const Cx cb = b.terms[0].second;
      for (const auto &[ma, ca] : a.terms) {
        const Cx c = ca * cb;
        if (!(c.re == 0 && c.im == 0)) r.terms.push_back({ma, c});
      }
      return r;
    }
    if (na == 1 && isConstant(a)) {
      NT_STAT_ADD(mul_const, 1);
      MPoly r(ns);
      r.terms.reserve(nb);
      const Cx ca = a.terms[0].second;
      for (const auto &[mb, cb] : b.terms) {
        const Cx c = ca * cb;
        if (!(c.re == 0 && c.im == 0)) r.terms.push_back({mb, c});
      }
      return r;
    }

    // Atom-free, inline-exponent operands (the Dirac fold and most contraction steps): a product
    // monomial is fully described by its packed exponent sum, so emit 24-byte keys plus a coefficient
    // array and build a Mono only for each surviving output term. The keys are emitted in the scratch
    // path's order and compared exactly as from_scratch compares atom-free inline monomials, so the
    // sort permutation and every like-term sum are identical. A field carry (an exponent > 31) needs
    // the heap representation and falls back to the scratch path below.
    if (na * nb <= kMulMaxScratch && MPoly::atomFreeInline(a) && MPoly::atomFreeInline(b)) {
      if (auto r = MPoly::mulKeyedAtomFree(a, b)) return std::move(*r);
    }
    // Operands carrying atoms: the same keyed product, with the atom multiset ranked (see mulKeyed).
    // Its setup only pays off on larger products; small ones stay on the scratch path.
    if (na * nb >= kMulKeyedMin && na * nb <= kMulMaxScratch) {
      if (auto r = MPoly::mulKeyed(a, b)) return std::move(*r);
    }

    if (na * nb <= kMulMaxScratch) { // exact byte-for-byte path — the common case
      MPolyScratch s;
      s.reserve(na * nb);
      for (const auto &[ma, ca] : a.terms)
        for (const auto &[mb, cb] : b.terms) emit(s, ma, ca, mb, cb);
      return MPoly::from_scratch(ns, std::move(s));
    }

    // Blocked: chunk the outer operand `a` so peak scratch is ~chunk·nb ≤ kMulMaxScratch, and fold the
    // per-chunk collapsed polynomials with operator+ (a linear merge that combines like terms again).
    NT_STAT_ADD(mul_blocked, 1);
    const std::size_t chunk = std::max<std::size_t>(1, kMulMaxScratch / nb);
    MPoly acc(ns);
    for (std::size_t i0 = 0; i0 < na; i0 += chunk) {
      const std::size_t i1 = std::min(na, i0 + chunk);
      MPolyScratch s;
      s.reserve((i1 - i0) * nb);
      for (std::size_t i = i0; i < i1; ++i) {
        const auto &[ma, ca] = a.terms[i];
        for (const auto &[mb, cb] : b.terms) emit(s, ma, ca, mb, cb);
      }
      MPoly part = MPoly::from_scratch(ns, std::move(s));
      acc = acc.terms.empty() ? std::move(part) : acc + part;
    }
    return acc;
  }

  inline MPoly MPoly::mulAdd(const MPoly &x0, const MPoly &y0, const MPoly &x1, const MPoly &y1)
  {
    // Only a sum of two non-trivial products gains: an empty or constant factor already has a
    // scratch-free path in operator*, and operator+ moves an empty side.
    auto keyable = [](const MPoly &x, const MPoly &y) {
      auto trivial = [](const MPoly &p) {
        if (p.terms.size() > 1) return false;
        if (p.terms.empty()) return true;
        const Mono &m = p.terms[0].first;
        return !m.e.overflow && m.e.packed[0] == 0 && m.e.packed[1] == 0 && m.atoms.empty();
      };
      return !trivial(x) && !trivial(y) && x.terms.size() * y.terms.size() <= kMulMaxScratch &&
             atomFreeInline(x) && atomFreeInline(y);
    };
    KeyedTerms t0, t1;
    if (!(keyable(x0, y0) && keyable(x1, y1) && keyedAtomFreeTerms(x0, y0, t0) && keyedAtomFreeTerms(x1, y1, t1)))
      return x0 * y0 + x1 * y1;
    const int ns = x0.nsym ? x0.nsym : (y0.nsym ? y0.nsym : (x1.nsym ? x1.nsym : y1.nsym));
    if (t0.empty()) return fromKeyed(ns, t1.data(), t1.size());
    if (t1.empty()) return fromKeyed(ns, t0.data(), t0.size());
    NT_STAT_ADD(add_calls, 1);
    KeyedTerms m; // operator+'s merge, in key form
    m.reserve(t0.size() + t1.size());
    std::size_t i = 0, j = 0;
    while (i < t0.size() && j < t1.size()) {
      const KeyedTerm &u = t0[i], &v = t1[j];
      if (u.p0 != v.p0 ? u.p0 < v.p0 : u.p1 < v.p1)
        m.push_back(t0[i++]);
      else if (v.p0 != u.p0 ? v.p0 < u.p0 : v.p1 < u.p1)
        m.push_back(t1[j++]);
      else {
        const Cx s = u.c + v.c;
        if (!(s.re == 0 && s.im == 0)) m.push_back({u.p0, u.p1, s});
        ++i;
        ++j;
      }
    }
    while (i < t0.size()) m.push_back(t0[i++]);
    while (j < t1.size()) m.push_back(t1[j++]);
    return fromKeyed(ns, m.data(), m.size());
  }

  /// @brief Cancel each numerator monomial against any atom whose denominator is a single monomial.
  ///
  /// `atomDen[aid]` is the atom's denominator polynomial `k²`. If it is a single monomial
  /// `c·∏ x^d`, then one power of the atom `1/D` cancels whenever the numerator exponents dominate
  /// `d` componentwise: subtract `d`, divide the coefficient by `c`, drop one atom instance. Repeat
  /// per atom occurrence. Atoms with a non-monomial denominator (shifted lines) are left intact here
  /// and handed to @ref divThroughPolyAtoms, which trial-divides them into the numerator; only an atom
  /// surviving both passes reaches the lowering as an `inv` env slot. Value-preserving, frame-agnostic.
  ///
  /// Split into a no-op TEST + a REBUILD so the `const MPoly&` and `MPoly&&` overloads below can share
  /// both and cannot drift apart. The `&&` overload exists because the pass-through case fires on ~95%
  /// of calls and `return p` on a `const&` is a full DEEP COPY of the polynomial — on the
  /// `eliminate` in-step reduction that is one copy per output index of every elimination step, which
  /// is the hottest loop in the engine. Same shape as the landed rvalue `operator+` lever.
  /// @brief The @ref divThroughMonomialAtoms pass-through test. True ⇒ the rebuild is the identity.
  inline bool dmaIsNoop(const MPoly &p, const std::vector<MPoly> &atomDen)
  {
    if (atomDen.empty()) return true;
    for (const auto &kv : p.terms)
      if (!kv.first.atoms.empty()) return false;
    return true; // no atom on any term (the common state after cancellation) ⇒ nothing can cancel
  }

  inline MPoly dmaRebuild(const MPoly &p, const std::vector<MPoly> &atomDen)
  {
    MPolyScratch out;
    out.reserve(p.terms.size());
    for (const auto &[m, c] : p.terms) {
      MonoExp e = m.e;
      Cx coeff = c;
      MonoAtoms keep;
      keep.reserve(m.atoms.size());
      for (const auto &aid : m.atoms) {
        bool cancelled = false;
        if (aid >= 0 && aid < (int)atomDen.size()) {
          const MPoly &D = atomDen[aid];
          if (D.terms.size() == 1) { // monomial denominator → may cancel
            const auto &dm = *D.terms.begin();
            const MonoExp &d = dm.first.e;
            const Cx dc = dm.second;
            bool dominates = dm.first.atoms.empty();
            for (int k = 0; dominates && k < p.nsym; ++k)
              if (e[k] < d[k]) dominates = false;
            if (dominates) {
              for (int k = 0; k < p.nsym; ++k)
                e[k] -= d[k];
              // coeff /= dc   (complex division: z / w = z·conj(w) / |w|²)
              const double den = dc.re * dc.re + dc.im * dc.im;
              assert(den != 0.0); // a stored monomial denominator never has a zero coefficient
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

  inline MPoly divThroughMonomialAtoms(const MPoly &p, const std::vector<MPoly> &atomDen)
  {
    NT_STAT_ADD(dma_calls, 1);
    if (dmaIsNoop(p, atomDen)) {
      NT_STAT_ADD(dma_noop, 1);
      return p;
    }
    return dmaRebuild(p, atomDen);
  }

  /// Rvalue overload: identical decision, but the pass-through case MOVES instead of deep-copying.
  inline MPoly divThroughMonomialAtoms(MPoly &&p, const std::vector<MPoly> &atomDen)
  {
    NT_STAT_ADD(dma_calls, 1);
    if (dmaIsNoop(p, atomDen)) {
      NT_STAT_ADD(dma_noop, 1);
      return std::move(p);
    }
    return dmaRebuild(p, atomDen);
  }

  /// @brief Apply unit-vector constraints `Σ_μ Uμ² = 1` to reduce each unit group's LAST component to
  ///        power ≤ 1: `U_last² → 1 − Σ_{μ<last} Uμ²`. A `group` is the list of symbol indices of one
  ///        unit direction's components (e.g. the loop direction `{U0,U1,U2}` with `U0²+U1²+U2²=1`, or a
  ///        `{cos,sin}` pair). This is the general form of `sin²→1−cos²` and is the crux of the numeric
  ///        backend's compactness: it (a) collapses the bare-loop denominator `Σ(l1·Uμ)² = l1²·ΣUμ²` to
  ///        the **monomial `l1²`** so @ref divThroughMonomialAtoms cancels it (like inv's `rel`), and
  ///        (b) collapses the `U·U` factors transverse projectors generate. Value-preserving;
  ///        terminates because each rewrite strictly lowers the last component's exponent.
  ///
  /// Split into a no-op TEST + a REBUILD so the `const MPoly&` and `MPoly&&` overloads below share
  /// both. The `&&` overload matters because the pass-through fires on ~88% of calls and `return p`
  /// on a `const&` deep-copies the whole polynomial; see @ref divThroughMonomialAtoms.
  /// @brief The @ref reduce_units pass-through test: no term carries any group's LAST component at
  ///        power >= 2, so the work-stack rebuild would be the identity (`p` is already sorted and
  ///        combined). Fires constantly — `eliminate`'s in-step reduction calls this on every
  ///        intermediate, most of which are already reduced.
  inline bool reduceUnitsIsNoop(const MPoly &p, const std::vector<std::vector<int>> &groups)
  {
    for (const auto &kv : p.terms)
      for (const auto &g : groups)
        if (!g.empty() && kv.first.e[g.back()] >= 2) return false;
    return true;
  }

  /// Keyed form of @ref reduceUnitsRebuild for inline exponents. The rewrite never touches a term's
  /// atoms, so a work item is its packed exponents plus the index of the source term whose atoms it
  /// carries, and a Mono is built only per output term. Pushes, pops, the output emission order, the
  /// sort comparator and the like-term sums are exactly those of the generic rebuild, so the result is
  /// bit-identical. Returns false (and leaves `r` unspecified) when an exponent would leave the inline
  /// range; the caller then runs the generic rebuild.
  inline bool reduceUnitsRebuildKeyed(const MPoly &p, const std::vector<std::vector<int>> &groups, MPoly &r)
  {
    for (const auto &g : groups)
      for (int k : g)
        if (k >= MonoExp::kInlineSyms) return false;
    for (const auto &kv : p.terms)
      if (kv.first.e.overflow) return false;
    auto field = [](const std::array<std::uint64_t, 2> &w, int k) {
      return static_cast<int>((w[MonoExp::wordOf(k)] >> MonoExp::shiftOf(k)) & MonoExp::kMaxExp);
    };
    auto unit = [](int k) { return std::uint64_t(1) << MonoExp::shiftOf(k); };
    struct Item {
      std::array<std::uint64_t, 2> e;
      std::uint32_t src;
      Cx c;
    };
    struct Key {
      std::uint64_t p0, p1;
      std::uint32_t idx;
    };
    std::vector<Item> work;
    work.reserve(p.terms.size());
    for (std::size_t t = 0; t < p.terms.size(); ++t)
      work.push_back({p.terms[t].first.e.packed, static_cast<std::uint32_t>(t), p.terms[t].second});
    std::vector<Key> keys;
    std::vector<std::uint32_t> outSrc;
    std::vector<Cx> outC;
    while (!work.empty()) {
      NT_STAT_ADD(ru_work, 1);
      const Item it = work.back();
      work.pop_back();
      int groupIdx = -1;
      for (int gIdx = 0; gIdx < (int)groups.size(); ++gIdx)
        if (!groups[gIdx].empty() && field(it.e, groups[gIdx].back()) >= 2) {
          groupIdx = gIdx;
          break;
        }
      if (groupIdx < 0) {
        keys.push_back({it.e[0], it.e[1], static_cast<std::uint32_t>(keys.size())});
        outSrc.push_back(it.src);
        outC.push_back(it.c);
        continue;
      }
      const std::vector<int> &group = groups[groupIdx];
      const int last = group.back();
      std::array<std::uint64_t, 2> base = it.e;
      base[MonoExp::wordOf(last)] -= 2 * unit(last);       // U_last^2 -> 1 - Σ_{μ<last} Uμ^2
      work.push_back({base, it.src, it.c});                // the "+1" branch
      for (std::size_t i = 0; i + 1 < group.size(); ++i) { // the "-Uμ^2" branches
        if (field(base, group[i]) + 2 > MonoExp::kMaxExp) return false;
        std::array<std::uint64_t, 2> shifted = base;
        shifted[MonoExp::wordOf(group[i])] += 2 * unit(group[i]);
        work.push_back({shifted, it.src, Cx{-it.c.re, -it.c.im}});
      }
    }
    auto atomsOf = [&](const Key &k) -> const MonoAtoms & { return p.terms[outSrc[k.idx]].first.atoms; };
    auto less = [&](const Key &x, const Key &y) {
      if (x.p0 != y.p0) return x.p0 < y.p0;
      if (x.p1 != y.p1) return x.p1 < y.p1;
      return atomsOf(x) < atomsOf(y);
    };
    bool sorted = true;
    for (std::size_t k = 1; sorted && k < keys.size(); ++k)
      sorted = less(keys[k - 1], keys[k]);
    if (!sorted) std::sort(keys.begin(), keys.end(), less);
    NT_STAT_ADD(fs_calls, 1);
    NT_STAT_ADD(fs_terms_in, keys.size());
    r = MPolyFactory::zero(p.nsym);
    r.terms.reserve(keys.size());
    for (const Key &k : keys) {
      const Cx c = outC[k.idx];
      if (c.re == 0 && c.im == 0) continue;
      if (!r.terms.empty()) {
        const Mono &back = r.terms.back().first;
        if (back.e.packed[0] == k.p0 && back.e.packed[1] == k.p1 && back.atoms == atomsOf(k)) {
          Cx &acc = r.terms.back().second;
          acc = acc + c;
          if (acc.re == 0 && acc.im == 0) r.terms.pop_back();
          continue;
        }
      }
      r.terms.emplace_back();
      r.terms.back().first.e.packed = {k.p0, k.p1};
      r.terms.back().first.atoms = atomsOf(k);
      r.terms.back().second = c;
    }
    return true;
  }

  inline MPoly reduceUnitsRebuild(const MPoly &p, const std::vector<std::vector<int>> &groups)
  {
    {
      MPoly r;
      if (reduceUnitsRebuildKeyed(p, groups, r)) return r;
    }
    MPolyScratch out;
    std::vector<std::tuple<MonoExp, MonoAtoms, Cx>> work;
    for (const auto &[m, c] : p.terms)
      work.push_back({m.e, m.atoms, c});
    while (!work.empty()) {
      NT_STAT_ADD(ru_work, 1);
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
      MonoExp base = e;
      base[last] -= 2;                                     // U_last^2 -> 1 - Σ_{μ<last} Uμ^2
      work.push_back({base, atoms, c});                    // the "+1" branch
      for (std::size_t i = 0; i + 1 < group.size(); ++i) { // the "-Uμ^2" branches
        MonoExp shifted = base;
        shifted[group[i]] += 2;
        work.push_back({std::move(shifted), atoms, Cx{-c.re, -c.im}});
      }
    }
    return MPoly::from_scratch(p.nsym, std::move(out));
  }

  /// Shared preamble of both @ref reduce_units overloads: the debug-only index validation. Kept as a
  /// separate function so the assert loop is written once; compiles away in release builds.
  inline void reduceUnitsCheckGroups([[maybe_unused]] const MPoly &p,
                                     [[maybe_unused]] const std::vector<std::vector<int>> &groups)
  {
    // every group entry is a symbol index, so it must address a valid component slot `e[idx]`
    for (const auto &g : groups)
      for ([[maybe_unused]] int idx : g)
        assert(idx >= 0 && idx < p.nsym);
  }

  inline MPoly reduce_units(const MPoly &p, const std::vector<std::vector<int>> &groups)
  {
    if (groups.empty()) return p; // NOT counted as a call, exactly as before
    NT_STAT_ADD(ru_calls, 1);
    reduceUnitsCheckGroups(p, groups);
    if (reduceUnitsIsNoop(p, groups)) {
      NT_STAT_ADD(ru_noop, 1);
      return p;
    }
    return reduceUnitsRebuild(p, groups);
  }

  /// Rvalue overload: identical decision, but both pass-through cases MOVE instead of deep-copying.
  inline MPoly reduce_units(MPoly &&p, const std::vector<std::vector<int>> &groups)
  {
    if (groups.empty()) return std::move(p); // NOT counted as a call, exactly as before
    NT_STAT_ADD(ru_calls, 1);
    reduceUnitsCheckGroups(p, groups);
    if (reduceUnitsIsNoop(p, groups)) {
      NT_STAT_ADD(ru_noop, 1);
      return std::move(p);
    }
    return reduceUnitsRebuild(p, groups);
  }

  /// @brief Cancel a MULTI-TERM inverse atom `1/D` against the numerator by exact polynomial division.
  ///
  /// @ref divThroughMonomialAtoms only cancels an atom whose denominator is a single **monomial**, and
  /// only term-by-term. A shifted line `k = l − q` has `k² = l² − 2 l·q + q²`, which is multi-term, so
  /// its atom always survived into an `inv` env slot — even when the numerator carries a factor of that
  /// very `k²`, which Dirac-trace numerators routinely do. This is the partial-fractioning step a
  /// fixed-frame contraction was assumed not to be able to do.
  ///
  /// Terms are grouped by their FULL atom multiset (so within a group every term carries the same
  /// atoms and the group is `G(x)·∏1/D_a`); each distinct multi-term `D_a` is then trial-divided into
  /// `G`. On an exact division `G·(1/D) → Q`, dropping one atom instance. Repeated to a fixed point.
  ///
  /// Exact and value-preserving: a division is only accepted when the remainder vanishes. "Vanishes"
  /// uses the same RELATIVE tolerance as the surrounding noise prune (`1e-9` against the dividend's
  /// largest coefficient) — a numeric frame makes exact cancellations land at ~1e-16 relative, far
  /// inside it, and the polynomial already carries round-off at that scale.
  ///
  /// Measured on ZAqbq1_147 Mq-in (real part): monomials 13,269 → 4,832, atom factors 1,168 → 832,
  /// fused SSA 33,775 → 7,649 (0.23x). 384 of 1,168 trial divisions are exact.
  /// Relative tolerance for "the division remainder vanishes". Same scale as the surrounding
  /// noise prune (`numeric_contract.hpp` `kNoisePruneRelTol`); kept local so `mpoly.hpp` stays
  /// independent of that header.
  inline constexpr double kPolyDivRelTol = 1e-9;

  /// The @ref divThroughPolyAtoms pass-through test: no term carries an atom whose denominator is
  /// multi-term with an atom-free leading term, so no trial division can run. The full pass then only
  /// regroups and re-sorts the (already canonical) terms and rebuilds each coefficient as `0 + c`.
  inline bool dpaIsNoop(const MPoly &p, const std::vector<MPoly> &atomDen)
  {
    for (const auto &kv : p.terms)
      for (const auto aid : kv.first.atoms)
        if (aid >= 0 && aid < (int)atomDen.size()) {
          const MPoly &D = atomDen[(std::size_t)aid];
          if (D.terms.size() >= 2 && D.terms.back().first.atoms.empty()) return false;
        }
    return true;
  }

  /// Rvalue overload: on a pass-through, normalise the coefficients in place (`0 + c` turns a −0.0
  /// component into +0.0, exactly as the full pass does) instead of regrouping and re-sorting.
  inline MPoly divThroughPolyAtoms(MPoly &&p, const std::vector<MPoly> &atomDen)
  {
    if (!dpaIsNoop(p, atomDen)) return divThroughPolyAtoms(static_cast<const MPoly &>(p), atomDen);
    NT_STAT_ADD(dpa_calls, 1);
    NT_STAT_ADD(dpa_noop, 1);
    for (auto &kv : p.terms) kv.second = Cx{0, 0} + kv.second;
    return std::move(p);
  }

  inline MPoly divThroughPolyAtoms(const MPoly &p, const std::vector<MPoly> &atomDen)
  {
    NT_STAT_ADD(dpa_calls, 1);
    // STEP 1 — group the polynomial by ATOM MULTISET. Only terms carrying the same set of inverse
    // atoms can cancel against the same denominator, and within one group the atoms are a common
    // factor that plays no part in the division. So each group is an ordinary polynomial in the
    // user symbols (`AtomGroup`), keyed by its exponent vector.
    using AtomGroup = std::map<MonoExp, Cx>;
    std::map<MonoAtoms, AtomGroup> byAtoms;
    for (const auto &[m, c] : p.terms) // monomials are unique, so each (atoms, e) slot is written once
      byAtoms[m.atoms].emplace(m.e, Cx{0, 0} + c);

    auto cdiv = [](const Cx &z, const Cx &w) { // complex divide; the coefficients are Cx, not double
      const double d = w.re * w.re + w.im * w.im;
      return Cx{(z.re * w.re + z.im * w.im) / d, (z.im * w.re - z.re * w.im) / d};
    };
    // Largest |coefficient| in a group — the scale the exactness test below is relative to. A trial
    // division is accepted only when its remainder is negligible against THIS, never against the
    // individual terms: two roundings of zero always disagree by 100%.
    auto maxAbsCoeff = [](const AtomGroup &G) {
      double m = 0.0;
      for (const auto &kv : G) m = std::max(m, std::max(std::fabs(kv.second.re), std::fabs(kv.second.im)));
      return m;
    };
    // Exact division G / D in lex order (the map's own ordering, a valid monomial order). Returns
    // false the moment a leading term is not divisible — i.e. the remainder is provably non-zero.
    auto divides = [&](AtomGroup G, const MPoly &D, AtomGroup &Q, double tol) {
      if (D.terms.size() < 2) return false;
      const auto &dlead = D.terms.back(); // MPoly::t is sorted by Mono, so back() is the lex-largest
      if (!dlead.first.atoms.empty()) return false;
      while (!G.empty()) {
        auto lead = std::prev(G.end());
        if (std::max(std::fabs(lead->second.re), std::fabs(lead->second.im)) < tol) { G.erase(lead); continue; }
        MonoExp e = lead->first;
        for (int k = 0; k < p.nsym; ++k) {
          if (e[k] < dlead.first.e[k]) return false;
          e[k] -= dlead.first.e[k];
        }
        const Cx c = cdiv(lead->second, dlead.second);
        Q[e] = Q[e] + c;
        for (const auto &dt : D.terms) {
          MonoExp f = e;
          for (int k = 0; k < p.nsym; ++k) f[k] += dt.first.e[k];
          Cx &g = G[f];
          g = Cx{g.re - (c.re * dt.second.re - c.im * dt.second.im),
                 g.im - (c.re * dt.second.im + c.im * dt.second.re)};
          if (std::max(std::fabs(g.re), std::fabs(g.im)) < tol) G.erase(f);
        }
      }
      return true;
    };

    MPolyScratch out;
    for (auto &[atoms0, G0] : byAtoms) {
      MonoAtoms atoms = atoms0;
      AtomGroup G = G0;
      for (bool again = true; again;) {
        again = false;
        const double tol = kPolyDivRelTol * maxAbsCoeff(G);
        for (std::size_t ai = 0; ai < atoms.size(); ++ai) {
          const int aid = atoms[ai];
          if (aid < 0 || aid >= (int)atomDen.size()) continue;
          const MPoly &D = atomDen[(std::size_t)aid];
          // Lead pre-filter: replicate EXACTLY the first check `divides` would make — its
          // multi-term/atom-free guards, then divisibility of the first significant (above-tol)
          // lead of G by D's lead — without paying the by-value map copy `divides` takes. A trial
          // failing on its first leading monomial is the common case (measured: >= 2/3 of trials,
          // and on some flows 100%). Trials that pass here still run `divides` unchanged, so the
          // accepted-division set and the resulting values are identical.
          if (D.terms.size() < 2 || !D.terms.back().first.atoms.empty()) {
            NT_STAT_ADD(dpa_pref, 1);
            continue;
          }
          {
            bool feasible = true;
            for (auto it = G.rbegin(); it != G.rend(); ++it) {
              if (std::max(std::fabs(it->second.re), std::fabs(it->second.im)) < tol) continue;
              const MonoExp &dl = D.terms.back().first.e;
              for (int k = 0; k < p.nsym; ++k)
                if (it->first[k] < dl[k]) {
                  feasible = false;
                  break;
                }
              break; // the first significant lead decides, exactly as divides' own loop would
            }
            if (!feasible) {
              NT_STAT_ADD(dpa_pref, 1);
              continue;
            }
          }
          AtomGroup Q;
          NT_STAT_ADD(dpa_trials, 1);
          if (divides(G, D, Q, tol)) {
            NT_STAT_ADD(dpa_exact, 1);
            G.swap(Q);
            atoms.erase(atoms.begin() + (long)ai);
            again = true;
            break;
          }
        }
      }
      for (auto &kv : G) {
        if (kv.second.re == 0.0 && kv.second.im == 0.0) continue;
        out.push_back({Mono{kv.first, atoms}, kv.second});
      }
    }
    return MPoly::from_scratch(p.nsym, std::move(out));
  }

  /// @brief Numeric evaluation (validation only). `x[i]` = user symbol i; `atomVal[aid]` = value of
  ///        `1/D_aid` (the caller supplies the reciprocal already evaluated).
  inline Cx eval(const MPoly &p, const std::vector<double> &x, const std::vector<double> &atomVal)
  {
    Cx s{0, 0};
    for (const auto &[m, c] : p.terms) {
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

  /// @brief Does every monomial of @p p carry an EVEN power of fundamental symbol @p sym?
  ///
  /// Used to certify a finite-T kernel as Matsubara-even, which lets DiFfRG's
  /// `QuadratureIntegrator_fT` replace `kernel(+ω) + kernel(−ω)` with `2·kernel(ω)` — halving both
  /// the emitted device code and the Matsubara-sum work at runtime.
  ///
  /// A monomial is `c · ∏ x_k^{e_k} · ∏ 1/atomDen_a`, so evenness of the whole polynomial in `sym`
  /// needs BOTH this predicate and the same predicate on every `atomDen` the polynomial references
  /// — an odd denominator would flip sign just as an odd numerator power does. The caller owns that
  /// second half (it is the one holding the atom table); this function deliberately does not guess
  /// at it.
  ///
  /// The condition is SUFFICIENT, not necessary: odd terms that cancel between monomials would be
  /// reported as odd. That is the safe direction — a false "even" silently drops the odd half of
  /// the Matsubara sum and gives wrong physics with no diagnostic, whereas a false "odd" only
  /// forgoes an optimisation.
  inline bool poly_even_in(const MPoly &p, int sym)
  {
    if (sym < 0) return false;
    for (const auto &[m, c] : p.terms) {
      (void)c;
      if (m.e.get(sym) % 2 != 0) return false;
    }
    return true;
  }

} // namespace numtracer::numeric
