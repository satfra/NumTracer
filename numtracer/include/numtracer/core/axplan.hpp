/// @file axplan.hpp
/// @brief The entry-type-agnostic per-axis contraction *planner* — the structural
///        index scheduler shared by the dense numeric tensor (@ref numtracer::dense::DTensor)
///        and any other labelled-axis tensor.
///
/// Each axis is a *type* @ref numtracer::core::Ax "Ax<Id, Dim>": `Id` is the
/// *contraction identity* (axes with the same `Id` on two tensors contract) and
/// `Dim` is that axis's extent. The planner is purely combinatorial — it resolves
/// which axes are shared (summed) versus free and the per-operand strides needed to
/// walk both operands — and carries no entry *values*, so it is independent of the
/// tensor's element type (expression types, `Cx`, ...). The element-specific fold
/// lives in the consuming tensor (the dense `contract`, ...).
#pragma once

#include "numtracer/core/config.hpp" // NT_THROW (exception-optional guard for -fno-exceptions builds)

#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace numtracer::core
{

  /// @brief Maximum tensor rank handled by the contraction plan scratch arrays.
  ///
  /// Must exceed the peak *intermediate* rank reached while a chain is left-folded,
  /// which can transiently outer-product (grow) before a shared axis closes it — not
  /// just the final rank. Larger networks (the 3-gluon Lorentz chains) need headroom
  /// here; an overflow silently corrupts the contraction.
  inline constexpr int kMaxAxisRank = 16;

  /// @brief Guard every @ref kMaxAxisRank-sized scratch write in the planners below.
  ///
  /// The doc on @ref kMaxAxisRank says an overflow "silently corrupts the contraction", and it used to: not
  /// one of the `[kMaxAxisRank]` arrays was bounds-checked. This is the *validation oracle*'s planner
  /// (@ref numtracer::dense::DTensor), so a silent corruption here does not merely produce a wrong
  /// number — it produces a wrong number that then CONFIRMS a wrong engine result. Loud is the only
  /// acceptable behaviour.
  ///
  /// Checking the operand ranks alone is NOT enough, which is the whole reason this exists as a
  /// named helper rather than two `if (RA > kMaxAxisRank)` lines: the result rank is `nFreeA + nFreeB`, so
  /// two perfectly in-range rank-10 operands **with no shared axis** produce `RR = 20 > kMaxAxisRank` and
  /// overrun `rid`/`rdim`. That is exactly the "transient outer-product growth" the @ref kMaxAxisRank doc
  /// warns about, and it is reachable from in-range inputs.
  ///
  /// `constexpr`-friendly: in a constant-evaluated call (@ref numtracer::dense::contract) a thrown
  /// exception is simply not a constant expression, so the overflow becomes a COMPILE error at the
  /// offending call site. At runtime (the `DynTensor` fold, `dense/dtensor.hpp`) it throws.
  constexpr void ax_check_rank(int r, const char *what)
  {
    if (r < 0 || r > kMaxAxisRank) NT_THROW(std::runtime_error, what);
  }

  /// @brief One tensor axis: a contraction identity and an extent.
  /// @tparam Id The contraction identity (axes with equal `Id` contract).
  /// @tparam Dim The axis's extent.
  template <int Id, int Dim> struct Ax {
    static constexpr int id = Id;   ///< Contraction identity.
    static constexpr int dim = Dim; ///< Axis extent.
  };

  /// @brief A compile-time list of axes; exposes per-axis id/dim arrays and total size.
  /// @tparam A The axis types (each an @ref Ax).
  template <class... A> struct AxList {
    static constexpr int rank = sizeof...(A);                                   ///< Number of axes.
    static constexpr std::array<int, (rank > 0 ? rank : 1)> ids = {A::id...};   ///< Per-axis identities.
    static constexpr std::array<int, (rank > 0 ? rank : 1)> dims = {A::dim...}; ///< Per-axis extents.
    static constexpr std::size_t size = (std::size_t{1} * ... * A::dim); ///< Product of dims.
  };

  // ---- mixed-radix helpers (per-axis dims) -----------------------------------

  /// @brief Decompose a flat index into per-axis indices (per-axis radices).
  /// @param flat The flat (row-major) index.
  /// @param dims The per-axis extents.
  /// @param rank The number of axes.
  /// @param out Output array; `out[a]` receives the index of axis `a`.
  constexpr void unflatten_mixed(std::size_t flat, const int *dims, std::size_t rank, std::size_t *out)
  {
    for (int a = rank - 1; a >= 0; --a) {
      out[a] = flat % dims[a];
      flat /= dims[a];
    }
  }
  /// @brief Compute row-major strides for per-axis extents.
  /// @param dims The per-axis extents.
  /// @param rank The number of axes.
  /// @param str Output array; `str[a]` receives the stride of axis `a`.
  constexpr void strides_mixed(const int *dims, int rank, std::size_t *str)
  {
    if (rank <= 0) return;
    str[rank - 1] = 1;
    for (int a = rank - 2; a >= 0; --a)
      str[a] = str[a + 1] * dims[a + 1];
  }

  // ---- contraction plan (structural => usable as a non-type template param) ----

  /// @brief A fully resolved contraction schedule for two labelled-axis tensors.
  ///
  /// A *structural type* so it can be passed as a non-type template parameter,
  /// which is how the result entries are built at compile time. It records which
  /// axes are shared (summed) versus free, the result axes, and the per-operand
  /// strides needed to walk both operands.
  struct EPlan {
    int RA = 0;                 ///< Rank of operand A.
    int RB = 0;                 ///< Rank of operand B.
    int RR = 0;                 ///< Result rank.
    int nSh = 0;                ///< Number of shared (summed) axes; entries occupy slots [0, nSh).
    int nFreeA = 0;             ///< Number of free axes of A.
    int nFreeB = 0;             ///< Number of free axes of B.
    int faAxis[kMaxAxisRank] = {};       ///< Free-axis positions of A.
    int fbAxis[kMaxAxisRank] = {};       ///< Free-axis positions of B.
    int saAxis[kMaxAxisRank] = {};       ///< Shared-axis positions in A.
    int sbAxis[kMaxAxisRank] = {};       ///< Shared-axis positions in B.
    int sdim[kMaxAxisRank] = {};         ///< Extent of each shared axis.
    int rid[kMaxAxisRank] = {};          ///< Identity of each result axis.
    int rdim[kMaxAxisRank] = {};         ///< Extent of each result axis.
    std::size_t strA[kMaxAxisRank] = {}; ///< Row-major strides of operand A.
    std::size_t strB[kMaxAxisRank] = {}; ///< Row-major strides of operand B.
    std::size_t nShared = 1;    ///< Product of the shared extents (the sum length).
    std::size_t nResult = 1;    ///< Product of the result extents (the entry count).
  };

  /// @brief Build the contraction plan for two axis lists.
  ///
  /// Pairs axes with matching identities as shared (summed) axes; the rest are free
  /// and become the result axes (A's first, then B's).
  /// @tparam Na Rank of operand A.
  /// @tparam Nb Rank of operand B.
  /// @param ida Operand A's per-axis identities.
  /// @param adim Operand A's per-axis extents.
  /// @param RA Rank of operand A.
  /// @param idb Operand B's per-axis identities.
  /// @param bdim Operand B's per-axis extents.
  /// @param RB Rank of operand B.
  /// @return The resolved @ref numtracer::core::EPlan.
  /// @warning Two axes that share an identity must agree on their extent; the plan
  ///          takes the shared dimension from A (assuming `adim == bdim` for the pair).
  template <std::size_t Na, std::size_t Nb>
  constexpr EPlan make_eplan(const std::array<int, Na> &ida, const std::array<int, Na> &adim, int RA,
                             const std::array<int, Nb> &idb, const std::array<int, Nb> &bdim, int RB)
  {
    EPlan p;
    // Guard BEFORE any [kMaxAxisRank] write. RA/RB bound the b_is_shared / faAxis / fbAxis / strA / strB
    // writes; p.RR (checked after it is known, below) bounds rid/rdim.
    ax_check_rank(RA, "make_eplan: operand A rank exceeds kMaxAxisRank (core/axplan.hpp)");
    ax_check_rank(RB, "make_eplan: operand B rank exceeds kMaxAxisRank (core/axplan.hpp)");
    p.RA = RA;
    p.RB = RB;
    bool b_is_shared[kMaxAxisRank] = {};
    for (int a = 0; a < RA; ++a) {
      int m = -1;
      for (int b = 0; b < RB; ++b)
        if (idb[b] == ida[a]) {
          m = b;
          break;
        }
      if (m >= 0) {
        // shared axes must agree on dimension — this was a COMMENT asserting an invariant nothing
        // checked (the `// == bdim[m]` below). A disagreeing pair silently contracts over A's
        // extent, reading past B or truncating it.
        if (adim[a] != bdim[m])
          NT_THROW(std::runtime_error,
                   "make_eplan: axes sharing an identity disagree on extent (core/axplan.hpp)");
        p.saAxis[p.nSh] = a;
        p.sbAxis[p.nSh] = m;
        p.sdim[p.nSh] = adim[a]; // == bdim[m]
        b_is_shared[m] = true;
        ++p.nSh;
      } else {
        p.faAxis[p.nFreeA++] = a;
      }
    }
    for (int b = 0; b < RB; ++b)
      if (!b_is_shared[b]) p.fbAxis[p.nFreeB++] = b;
    p.RR = p.nFreeA + p.nFreeB;
    // THE case rank-checking the operands misses: two in-range operands sharing NO axis
    // outer-product to RR = RA + RB. Two rank-10 operands give RR = 20 > kMaxAxisRank and overrun rid/rdim.
    ax_check_rank(p.RR, "make_eplan: result rank nFreeA+nFreeB exceeds kMaxAxisRank — the chain "
                        "outer-products before a shared axis closes it (core/axplan.hpp)");
    for (int t = 0; t < p.nFreeA; ++t) {
      p.rid[t] = ida[p.faAxis[t]];
      p.rdim[t] = adim[p.faAxis[t]];
    }
    for (int t = 0; t < p.nFreeB; ++t) {
      p.rid[p.nFreeA + t] = idb[p.fbAxis[t]];
      p.rdim[p.nFreeA + t] = bdim[p.fbAxis[t]];
    }
    // The arrays may be padded to size 1 for rank-0 operands; loops use RA/RB only.
    int a_dims_local[kMaxAxisRank] = {}, b_dims_local[kMaxAxisRank] = {};
    for (int a = 0; a < RA; ++a)
      a_dims_local[a] = adim[a];
    for (int b = 0; b < RB; ++b)
      b_dims_local[b] = bdim[b];
    strides_mixed(a_dims_local, RA, p.strA);
    strides_mixed(b_dims_local, RB, p.strB);
    p.nShared = 1;
    for (int k = 0; k < p.nSh; ++k)
      p.nShared *= p.sdim[k];
    p.nResult = 1;
    for (int i = 0; i < p.RR; ++i)
      p.nResult *= p.rdim[i];
    return p;
  }

  /// @brief Flat index into operand A for a given result index and shared index.
  /// @param r The flat result index.
  /// @param s The flat shared (summed) index.
  /// @param p The contraction plan.
  /// @return The flat index into A's entries.
  constexpr std::size_t a_index(std::size_t r, std::size_t s, const EPlan &p)
  {
    std::size_t cr[kMaxAxisRank] = {}, cs[kMaxAxisRank] = {};
    unflatten_mixed(r, p.rdim, p.RR, cr);
    unflatten_mixed(s, p.sdim, p.nSh, cs);
    std::size_t a = 0;
    for (int t = 0; t < p.nFreeA; ++t)
      a += cr[t] * p.strA[p.faAxis[t]];
    for (int k = 0; k < p.nSh; ++k)
      a += cs[k] * p.strA[p.saAxis[k]];
    return a;
  }
  /// @brief Flat index into operand B for a given result index and shared index.
  /// @param r The flat result index.
  /// @param s The flat shared (summed) index.
  /// @param p The contraction plan.
  /// @return The flat index into B's entries.
  constexpr std::size_t b_index(std::size_t r, std::size_t s, const EPlan &p)
  {
    std::size_t cr[kMaxAxisRank] = {}, cs[kMaxAxisRank] = {};
    unflatten_mixed(r, p.rdim, p.RR, cr);
    unflatten_mixed(s, p.sdim, p.nSh, cs);
    std::size_t b = 0;
    for (int t = 0; t < p.nFreeB; ++t)
      b += cr[p.nFreeA + t] * p.strB[p.fbAxis[t]];
    for (int k = 0; k < p.nSh; ++k)
      b += cs[k] * p.strB[p.sbAxis[k]];
    return b;
  }

  /// @brief Build the result @ref numtracer::core::AxList type from a plan's result axes.
  /// @tparam plan The contraction plan.
  /// @tparam I The result-axis index pack `0..RR-1`.
  /// @return (unevaluated) the `AxList<Ax<rid,rdim>...>` result-axis type.
  template <EPlan plan, std::size_t... I>
  auto make_axlist(std::index_sequence<I...>) -> AxList<Ax<plan.rid[I], plan.rdim[I]>...>;

  // ---- entrywise sum plan (the dual of contract) ------------------------------

  /// @brief A resolved plan for adding two tensors with the same axis-id set.
  ///
  /// Addition needs both operands over the *same* axes, but the builders may produce
  /// them in different axis orders. This records, for each axis position of A, the
  /// row-major stride of the B axis carrying the same identity, so B is reindexed into
  /// A's layout — the additive analogue of @ref numtracer::core::EPlan.
  struct EAddPlan {
    int R = 0;                  ///< Rank (shared by both operands and the result).
    int rdim[kMaxAxisRank] = {};         ///< Result (A's) per-axis extents.
    std::size_t bstr[kMaxAxisRank] = {}; ///< `bstr[a]` = stride of B's axis whose id equals A's axis `a`.
    std::size_t n = 1;          ///< Entry count (product of dims).
  };
  /// @brief Build the add plan: match axes by identity, take B's stride for each A axis.
  /// @warning A and B must carry the same set of axis identities with matching extents.
  template <std::size_t Na, std::size_t Nb>
  constexpr EAddPlan make_eaddplan(const std::array<int, Na> &ida, const std::array<int, Na> &adim, int R,
                                   const std::array<int, Nb> &idb, const std::array<int, Nb> &bdim)
  {
    EAddPlan p;
    ax_check_rank(R, "make_eaddplan: rank exceeds kMaxAxisRank (core/axplan.hpp)");
    p.R = R;
    int b_dims_local[kMaxAxisRank] = {};
    std::size_t bstrideAll[kMaxAxisRank] = {};
    for (int b = 0; b < R; ++b)
      b_dims_local[b] = bdim[b];
    strides_mixed(b_dims_local, R, bstrideAll);
    for (int a = 0; a < R; ++a) {
      p.rdim[a] = adim[a];
      for (int b = 0; b < R; ++b)
        if (idb[b] == ida[a]) {
          p.bstr[a] = bstrideAll[b];
          break;
        } // B axis with the same id
      p.n *= adim[a];
    }
    return p;
  }
  /// @brief B's flat index for result index `r` (which is in A's layout).
  constexpr std::size_t add_bindex(std::size_t r, const EAddPlan &p)
  {
    std::size_t cr[kMaxAxisRank] = {};
    unflatten_mixed(r, p.rdim, p.R, cr);
    std::size_t b = 0;
    for (int a = 0; a < p.R; ++a)
      b += cr[a] * p.bstr[a];
    return b;
  }

} // namespace numtracer::core
