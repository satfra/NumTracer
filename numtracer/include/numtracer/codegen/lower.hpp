/// @file lower.hpp
/// @brief The Horner core that lowers a set of monomials (@ref LMono) to a real straight-line
///        program through the shared CSE back-end.
///
/// A contracted diagram polynomial is a set of monomials, each a real coefficient times a product
/// of `var(envId)^exp` factors (the env id names a fundamental symbol — a scalar product, an
/// inverse propagator, a dressing call, or a raw kernel argument). Evaluating that fast is the job
/// of hash-consed value-numbering (@ref numtracer::network::rdetail::RBuilder — the CSE) plus greedy
/// multivariate Horner factoring: @ref horner emits the monomial set into an `RBuilder`, sharing
/// powers / products / factors, and returns the result slot. The numeric backend
/// (@ref numtracer::numeric::to_genprog) builds the @ref LMono set and drives this via
/// @ref numtracer::network::gdetail::best_into in `codegen/gen.hpp`.
#pragma once

#include "numtracer/codegen/real_cse.hpp"

#include <algorithm>
#include <climits>
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

/// @brief Greedy multivariate Horner of a monomial set, emitted through the real value-numbering
///        builder (so shared powers / products / factors are CSE'd). Returns the result slot.
constexpr int horner(rdetail::RBuilder &builder, std::vector<LMono> terms) {
  using namespace rdetail;
  // ---- base case: all monomials are constant → emit their sum -----------------------------
  bool allconst = true;
  for (const LMono &m : terms)
    if (!m.vp.empty()) { allconst = false; break; }
  if (allconst) {
    double s = 0;
    for (const LMono &m : terms) s += m.c;
    return rconst(builder, s);
  }
  const auto [pivot, pivotExp] = choose_pivot(terms);
  auto [with, without] = partition_pivot(std::move(terms), pivot, pivotExp);
  // ---- recurse on both parts, then combine: pivot^pivotExp * horner(with) + horner(without) -----
  const int withSlot = horner(builder, std::move(with));
  const int withoutSlot = without.empty() ? -1 : horner(builder, std::move(without));
  const int pivotVar = rvar(builder, pivot);
  int pivotPow = pivotVar;
  for (int e = 1; e < pivotExp; ++e) pivotPow = rmul(builder, pivotPow, pivotVar);
  const int factored = rmul(builder, pivotPow, withSlot);
  return radd(builder, factored, withoutSlot);
}

} // namespace numtracer::network
