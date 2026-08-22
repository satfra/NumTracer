/// @file lower.hpp
/// @brief The Horner core that lowers a set of monomials (@ref LMono) to a real straight-line
///        program through the shared CSE back-end.
///
/// A contracted diagram polynomial is a set of monomials, each a real coefficient times a product
/// of `var(envId)^exp` factors (the env id names a fundamental symbol — a scalar product, an
/// inverse propagator, a dressing call, or a raw kernel argument). Evaluating that fast is the job
/// of hash-consed value-numbering (@ref numtracer::network::rdetail::RBuilder — the CSE) plus greedy
/// multivariate Horner factoring: @ref horner emits the monomial set into an `RBuilder`, sharing
/// powers / products / factors.
///
/// It returns an @ref NVal, not a bare slot: the monomial COEFFICIENT is deliberately kept out of the
/// emitted graph and handed back to the caller, because the CSE keys constants on their raw bits and a
/// shape reached under a different coefficient would otherwise intern separately and be recomputed.
/// @ref numtracer::network::gdetail::best_into in `codegen/gen.hpp` is what materialises that scalar
/// (once per trace) and hands an ordinary slot to the numeric backend
/// (@ref numtracer::numeric::to_genprog), which is also what builds the @ref LMono set.
#pragma once

#include "numtracer/codegen/real_cse.hpp"
#include "numtracer/core/envvar.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <utility>
#include <vector>

namespace numtracer::network {

/// @brief A monomial of the lowered polynomial: `c * prod(var(envId)^exp)`, `vp` sorted by envId.
struct LMono {
  double c = 0;
  std::vector<std::pair<int, int>> vp; ///< (envId, exponent)
};

/// @brief Pick the Horner pivot: the variable occurring in the most terms, and the lowest exponent
///        `pivotExp` it appears with (the most we can factor out). Returns `{pivot envId, pivotExp}`.
///
///        The tally is a DENSE id-indexed count (env ids are small contiguous ints) plus a
///        first-appearance-ordered id list: O(Σ|vp| + V) instead of the old per-factor linear scan's
///        O(Σ|vp|·V), which at V≈40 distinct ids was the dominant term of the whole lowering
///        (O(N·|vp|·V²) per trace). Ties MUST resolve to the id that first appears walking the terms
///        in order — that reproduces the old insertion-order scan exactly, and any other tie-break
///        (e.g. smallest id) changes pivots and thus the emitted kernel bytes.
constexpr std::pair<int, int> choose_pivot(const std::vector<LMono> &terms) {
  std::vector<int> count;  // indexed by envId
  std::vector<int> order;  // envIds in first-appearance order (the tie-break)
  for (const LMono &m : terms)
    for (auto [id, e] : m.vp) {
      if (id >= (int)count.size()) count.resize((std::size_t)id + 1, 0);
      if (count[(std::size_t)id]++ == 0) order.push_back(id);
    }
  int pivot = -1, bestCount = -1;
  for (int id : order)
    if (count[(std::size_t)id] > bestCount) { bestCount = count[(std::size_t)id]; pivot = id; }
  int pivotExp = INT_MAX;
  for (const LMono &m : terms)
    for (auto [id, ex] : m.vp)
      if (id == pivot) pivotExp = std::min(pivotExp, ex);
  return {pivot, pivotExp};
}

/// @brief Partition @p terms into `with` (those containing `pivot`, with `pivot^pivotExp` divided
///        out) and `without` (the rest). Returns `{with, without}`.
constexpr std::pair<std::vector<LMono>, std::vector<LMono>>
partition_pivot(std::vector<LMono> terms, int pivot, int pivotExp) {
  std::vector<LMono> with, without;
  for (LMono &m : terms) {
    // Divide `pivot^pivotExp` out IN PLACE — adjust or erase the pivot entry and move the monomial,
    // reusing its vp buffer, instead of rebuilding a fresh vector per monomial per recursion level
    // (that was O(N·V) heap allocations per trace). Entry order is untouched by the in-place edit
    // (entries stay sorted by envId), so the result — and the emitted kernel — is byte-identical.
    bool has = false;
    for (std::size_t i = 0; i < m.vp.size(); ++i)
      if (m.vp[i].first == pivot) {
        if (m.vp[i].second > pivotExp) m.vp[i].second -= pivotExp;
        else m.vp.erase(m.vp.begin() + (std::ptrdiff_t)i);
        has = true;
        break;
      }
    if (has) with.push_back(std::move(m));
    else without.push_back(std::move(m));
  }
  return {std::move(with), std::move(without)};
}

/// @brief A lowered value in *normalised* form: `value == g * (s < 0 ? 1.0 : slot s)`.
///
/// The point of carrying the scalar separately is that the SSA value numbering
/// (@ref numtracer::network::rdetail::RBuilder) keys constants on their raw bits, so the SAME
/// expression shape reached with a different coefficient interns as a different node and is
/// recomputed. Emitting `295.6*(f31 + f29*f30)` with the constant pushed into the leaves — which is
/// what a plain Horner does — makes `(f31 + f29*f30)` unshareable with the dozens of other places
/// that same shape occurs under a different coefficient. Keeping the scalar out of the graph and
/// folding it into the PARENT's ratio instead lets the shape intern once.
///
/// `s < 0` means "the value is exactly the scalar `g`", which subsumes the structural zero
/// (`g == 0`), the structural one, and any pure-constant subtree — so the all-constant base case
/// emits no instruction at all.
///
/// DANGER: `rdetail::rmul` / `rdetail::radd` read a negative slot as structural ZERO. An `NVal`
/// with `s < 0` therefore must never be handed to them directly; go through @ref scale_into or
/// @ref mul_shape, which is why neither the type nor a raw `.s` escapes this header.
struct NVal {
  double g = 0;   ///< the scalar.
  int s = -1;     ///< slot of the constant-free shape, or < 0 when the value is just `g`.
};

/// @brief Materialise an @ref NVal as an ordinary slot (structural zero stays the `-1` sentinel).
inline int scale_into(rdetail::RBuilder &builder, NVal v)
{
  using namespace rdetail;
  if (v.s < 0) return rconst(builder, v.g); // covers g == 0 -> -1
  if (v.g == 1.0) return v.s;
  return rmul(builder, rconst(builder, v.g), v.s); // rmul folds x1 and x(-1) itself
}

/// @brief Wrap an already-materialised slot as an @ref NVal, mapping the structural-zero sentinel
///        onto the `{0.0, -1}` form the invariant requires.
///
/// Needed because `-1` means two different things on the two sides of this type: to the builder it is
/// structural ZERO, to an `NVal` shape it is the structural ONE. `{1.0, -1}` would therefore read back
/// as the constant 1 -- which is exactly how a polynomial that cancelled to nothing turned into a
/// kernel returning 1.0 instead of 0.0.
inline NVal from_slot(int s) { return s < 0 ? NVal{0.0, -1} : NVal{1.0, s}; }

/// @brief `pivotPow * shape`, where a shape of `s < 0` is the structural one.
inline int mul_shape(rdetail::RBuilder &builder, int pivotPow, int shape)
{
  return shape < 0 ? pivotPow : rdetail::rmul(builder, pivotPow, shape);
}

/// @brief Is the scalar-normalising lowering enabled? `NT_GEN_NO_NORMHORNER` turns it off.
inline bool normhorner_enabled()
{
  static const bool on = !env_flag("NT_GEN_NO_NORMHORNER");
  return on;
}

/// @brief Greedy multivariate Horner of a monomial set, emitted through the real value-numbering
///        builder (so shared powers / products / factors are CSE'd).
///
/// With @p normalise the result is returned in @ref NVal form — the coefficient is kept OUT of the
/// emitted graph and handed to the caller, which folds it into its own ratio. With it false this
/// reproduces the classic lowering exactly, returning `{1.0, slot}`; `best_into` costs both and keeps
/// the smaller, because on a trace with no repeated shapes normalising buys nothing and still
/// restructures.
///
/// A common-factor extraction at the combine (`A*X + A*Y*Z -> A*(X+Y*Z)`) was implemented and
/// REJECTED: at a combine the two addends are `pivot^n * W` and `ratio * Wo`, whose only factors are
/// the pivot chain (absent from the other by construction) and the branch sums themselves, so a shared
/// factor needs `W == Wo`. Measured 1.000x/0.996x/0.999x ops on ZAqbq{1,4,7}_147 and no runtime gain on
/// either CPU or GPU. The structure it wants is real but lives elsewhere — 27% of add nodes on the
/// production with_mesons/ZA4 kernel versus 4-9% here — so it would need a post-pass over the finished
/// instruction stream, not a hook at this combine.
inline NVal horner(rdetail::RBuilder &builder, std::vector<LMono> terms, bool normalise)
{
  using namespace rdetail;
  // ---- base case: all monomials are constant → their sum, with no instruction emitted ----------
  bool allconst = true;
  for (const LMono &m : terms)
    if (!m.vp.empty()) { allconst = false; break; }
  if (allconst) {
    double sum = 0;
    for (const LMono &m : terms) sum += m.c;
    // The un-normalised arm keeps materialising the constant HERE, at the base, rather than deferring
    // it to the parent's combine: that is what makes it emit byte-identically to the classic lowering,
    // so `NT_GEN_NO_NORMHORNER` is a true A/B control and not merely an equivalent one.
    if (!normalise) return from_slot(rconst(builder, sum));
    return {sum, -1};
  }
  const auto [pivot, pivotExp] = choose_pivot(terms);
  auto [with, without] = partition_pivot(std::move(terms), pivot, pivotExp);
  // ---- recurse on both parts, then combine: pivot^pivotExp * horner(with) + horner(without) -----
  const NVal W = horner(builder, std::move(with), normalise);
  const NVal Wo = without.empty() ? NVal{0.0, -1} : horner(builder, std::move(without), normalise);
  const int pivotVar = rvar(builder, pivot);
  int pivotPow = pivotVar;
  for (int e = 1; e < pivotExp; ++e) pivotPow = rmul(builder, pivotPow, pivotVar);

  if (!normalise)
    return from_slot(radd(builder, rmul(builder, pivotPow, scale_into(builder, W)),
                          scale_into(builder, Wo)));

  if (W.g == 0.0) return Wo;                       // the with-branch cancelled away
  const int P = mul_shape(builder, pivotPow, W.s); // the with-branch shape, coefficient stripped
  if (Wo.g == 0.0) return {W.g, P};

  // Keep the LARGER-magnitude scalar, so the emitted ratio has |ratio| <= 1. That choice is
  // scale-covariant — scaling a subtree's coefficients by k scales both branch scalars, so the same
  // branch still wins and the emitted ratio is unchanged — hence the CSE key is exactly as canonical
  // as always keeping the with-branch's scalar would be, and it costs no instruction (the fma
  // absorbs one multiply either way). It also makes a magnitude window unnecessary: such a window
  // would fall back to the un-normalised form for individual nodes, whose shapes then fail to match
  // their normalised siblings, fragmenting the sharing exactly where coefficients spread widest.
  //
  // Exception: when the without-branch is a bare constant (`Wo.s < 0`) its shape is the structural
  // one, and keeping ITS scalar would force a literal `1` into the graph. Keep the with-branch's
  // scalar there regardless of magnitude.
  const bool keepW = Wo.s < 0 || std::fabs(W.g) >= std::fabs(Wo.g);
  const double g = keepW ? W.g : Wo.g;
  // NOT rounded. Rounding the ratio would merge coefficient noise (`-31.999999999999002` -> `-32`)
  // and buys ~10% more sharing, but the noise is ~1e-14 relative, so any rounding coarse enough to
  // merge it perturbs by ~1e-14 — and unlike `snap_coeff`, which acts once per monomial at a leaf,
  // that lands at EVERY interior node where the ~1e6 cancellation amplification acts. Measured on
  // ZA3_147 against compare_za3_147_num's 1e-9 gate: unrounded 5.6e-11 (healthy), rounded to 13
  // significant digits 2.7e-08 (FAILING). No intermediate setting helps, because the noise and the
  // harmful perturbation are the same size; recovering that 10% needs exact rational coefficients
  // upstream, not a rounding here.
  const double ratio = (keepW ? Wo.g : W.g) / g;
  // A ratio that underflowed to zero would delete the other branch outright (rconst(0) is the
  // structural-zero sentinel, which rmul propagates and radd then drops); one that overflowed would
  // emit an infinity. Neither is reachable once |ratio| <= 1 and the relative noise prune bounds the
  // coefficient spread, but both are silent-wrong-answer failures, so check rather than assume.
  if (!std::isfinite(ratio) || ratio == 0.0)
    return {1.0, radd(builder, rmul(builder, rconst(builder, W.g), P), scale_into(builder, Wo))};

  const int rc = rconst(builder, ratio);
  if (keepW) {
    // value == g * (P + ratio * Wo_shape); a constant without-branch degenerates to the ratio itself
    const int addend = Wo.s < 0 ? rc : rmul(builder, rc, Wo.s);
    return {g, radd(builder, P, addend)};
  }
  return {g, radd(builder, rmul(builder, rc, P), Wo.s)};
}

} // namespace numtracer::network
