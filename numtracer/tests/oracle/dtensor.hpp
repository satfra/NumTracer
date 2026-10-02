/// @file dtensor.hpp
/// @brief Test-side **dense** numeric tensor: a brute-force oracle that contracts labelled axes
///        entry-for-entry. Used by test_dense_pack, test_guards and the two frozen dense benchmarks
///        (bench_za_dense, bench_za3_147_dense); not part of the installed library.
///
/// `DTensor<AX>` carries per-axis @ref numtracer::core::Ax "Ax<Id,Dim>" labels (so one tensor
/// can mix Lorentz(4) ⊗ Dirac-spinor(4) ⊗ SU(N) at once), with its entries runtime
/// @ref numtracer::Cx **numbers** in a `std::array`. The type depends only on the `AxList`, so a
/// dense kernel compiles flat, and contraction sums **every** index combination (no
/// structural-zero pruning). The contraction schedule comes from the planner in `axplan.hpp`.
///
/// Builders: Lorentz `metric`/`vector`/`transverse_projector`/`longitudinal_projector`, Dirac
/// `gamma_axis`/`identity`, colour `T`/`f`/`delta_fund`/`delta_adj`, filled from the library's
/// data tables and the `double renv[]` runtime environment.
#pragma once

#include "axplan.hpp"                     // Ax, AxList, EPlan + planner
#include "numtracer/core/cx.hpp"          // Cx
#include "numtracer/dirac/dirac_data.hpp" // gamma_entry
#include "numtracer/sun/sun_data.hpp"     // SUNData<N> (generators, f_nonzeros)

#include <array>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

// Region guard disabling loop auto-vectorization: GCC 16.1.1 -O3 -march=native (AVX-512VL) miscompiles
// the complex-`Cx` std::vector fold loops below (empty result buffer -> SIGSEGV). Clang is unaffected.
#if defined(__GNUC__) && !defined(__clang__)
#define NT_BEGIN_NO_LOOP_VECTORIZE _Pragma("GCC push_options") _Pragma("GCC optimize(\"no-tree-loop-vectorize\")")
#define NT_END_NO_LOOP_VECTORIZE _Pragma("GCC pop_options")
#else
#define NT_BEGIN_NO_LOOP_VECTORIZE
#define NT_END_NO_LOOP_VECTORIZE
#endif

namespace numtracer::dense
{

  using core::Ax;
  using core::AxList;

  /// @brief A dense numeric tensor over labelled axes — the entry-for-entry numeric oracle.
  /// @tparam AX The axis list (an @ref numtracer::core::AxList).
  template <class AX> struct DTensor {
    using axes = AX;                 ///< The axis list.
    std::array<Cx, AX::size> data{}; ///< Row-major numeric entries (all of them).
  };

  // ---- allocation-lean dynamic-rank fold (the variadic contract_all/add_all/scale path) -------
  //
  // `DTensor<AX>` stores its entries inline in `std::array<Cx, AX::size>`. In a long contraction
  // chain the running intermediate can transiently outer-product to a high rank, so that array
  // grows to megabytes — and the whole chain is materialised in one inlined stack frame, enough
  // to overflow the stack on the largest gluon traces (a single ZA3 1/4/7 trace builds a multi-GB
  // frame). The fold below instead keeps the running intermediate on the heap (a `DynTensor`'s
  // `std::vector<Cx>`) so any rank is just a small handle on the stack, and it is careful about
  // heap traffic in two ways:
  //   - operands are *viewed in place* (@ref detail::View) — a typed builder is contracted
  //     straight out of its inline `std::array`, never copied onto the heap;
  //   - the running intermediate *ping-pongs between two reusable buffers* (`ping`/`pong`), so a
  //     chain of N factors performs O(1) allocations (the two buffers, one of which is moved out
  //     as the result) instead of one per step — and the big transient tensors are reused, not
  //     reallocated.
  // The arithmetic is the identical brute-force contraction (same planner, same serial summation
  // order, same numbers) — only the storage discipline changes. All members/functions are
  // `constexpr`, so the small colour-only "constant" traces some kernels evaluate at compile time
  // keep folding (the transient vector allocations are reclaimed within the constant expression).

  /// @brief A heap-backed dynamic-rank numeric tensor — the result carrier of the variadic fold,
  ///        and the reusable scratch buffer it ping-pongs through. `data[0]` is the scalar result
  ///        of a closed loop.
  struct DynTensor {
    std::vector<Cx> data;              ///< Heap-backed row-major entries (all of them).
    int rank = 0;                      ///< Number of axes.
    std::array<int, core::kMaxAxisRank> ids{};  ///< Per-axis identities (first @ref rank valid).
    std::array<int, core::kMaxAxisRank> dims{}; ///< Per-axis extents (first @ref rank valid).
    /// @brief The rank-0 (scalar) result as a `Cx`.
    constexpr Cx scalar_value() const { return data[0]; }
    /// @brief The rank-0 (scalar) result as a `std::complex<double>`.
    constexpr std::complex<double> scalar_cplx() const { return {data[0].re, data[0].im}; }
  };

  namespace detail
  {
    /// @brief A non-owning view of a tensor's entries + axis metadata. Lets the fold treat a typed
    ///        builder @ref DTensor and a dynamic intermediate @ref DynTensor uniformly, without
    ///        copying either onto the heap.
    struct View {
      const Cx *data = nullptr;  ///< Row-major entries (borrowed).
      int rank = 0;              ///< Number of axes.
      const int *ids = nullptr;  ///< Per-axis identities (borrowed, @ref rank valid).
      const int *dims = nullptr; ///< Per-axis extents (borrowed, @ref rank valid).
    };
  } // namespace detail

  /// @brief View a typed builder in place (its entries live in the inline `std::array`).
  template <class AX> constexpr detail::View view_of(const DTensor<AX> &t)
  {
    return {t.data.data(), AX::rank, AX::ids.data(), AX::dims.data()};
  }
  /// @brief View a dynamic intermediate in place.
  constexpr detail::View view_of(const DynTensor &d) { return {d.data.data(), d.rank, d.ids.data(), d.dims.data()}; }

// No loop vectorization for the fold block (the GCC miscompile noted at NT_BEGIN_NO_LOOP_VECTORIZE).
NT_BEGIN_NO_LOOP_VECTORIZE

  /// @brief Contract views `a`,`b` over their shared axis identities, writing the result into
  ///        `out` (its `data` vector is resized, reusing capacity across fold steps). `out` must
  ///        not alias either operand's storage.
  constexpr void contract_into(DynTensor &out, detail::View a, detail::View b)
  {
    std::array<int, core::kMaxAxisRank> a_ids{}, a_dims{}, b_ids{}, b_dims{};
    for (int i = 0; i < a.rank; ++i) {
      a_ids[i] = a.ids[i];
      a_dims[i] = a.dims[i];
    }
    for (int i = 0; i < b.rank; ++i) {
      b_ids[i] = b.ids[i];
      b_dims[i] = b.dims[i];
    }
    // Build the contraction plan (also yields the result tensor's layout).
    const core::EPlan plan = core::make_eplan(a_ids, a_dims, a.rank, b_ids, b_dims, b.rank);
    out.rank = plan.RR;
    for (int t = 0; t < plan.RR; ++t) {
      out.ids[t] = plan.rid[t];
      out.dims[t] = plan.rdim[t];
    }
    out.data.resize(plan.nResult);
    // Division-free odometers: instead of recomputing a mixed-radix decomposition (a modulo/divide
    // per axis) of r and s at every step — which dominated the runtime — precompute each axis's
    // flat-offset increment and walk the result- and shared-index counters by add/sub with carry.
    std::size_t result_stride_a[core::kMaxAxisRank] = {}, result_stride_b[core::kMaxAxisRank] = {}; // result axis t -> stride into A / B (one is 0)
    for (int t = 0; t < plan.nFreeA; ++t)
      result_stride_a[t] = plan.strA[plan.faAxis[t]];
    for (int t = 0; t < plan.nFreeB; ++t)
      result_stride_b[plan.nFreeA + t] = plan.strB[plan.fbAxis[t]];
    std::size_t shared_stride_a[core::kMaxAxisRank] = {}, shared_stride_b[core::kMaxAxisRank] = {}; // shared axis k -> stride into A / B
    for (int k = 0; k < plan.nSh; ++k) {
      shared_stride_a[k] = plan.strA[plan.saAxis[k]];
      shared_stride_b[k] = plan.strB[plan.sbAxis[k]];
    }

    int cur_result_idx[core::kMaxAxisRank] = {};
    std::size_t baseA = 0, baseB = 0;
    for (std::size_t r = 0; r < plan.nResult; ++r) {
      int cur_shared_idx[core::kMaxAxisRank] = {};
      std::size_t offA = baseA, offB = baseB;
      Cx acc{0, 0};
      for (std::size_t s = 0; s < plan.nShared; ++s) {
        acc = acc + a.data[offA] * b.data[offB];
        for (int k = plan.nSh - 1; k >= 0; --k) { // step the shared odometer (no-op when nSh == 0)
          if (++cur_shared_idx[k] < plan.sdim[k]) {
            offA += shared_stride_a[k];
            offB += shared_stride_b[k];
            break;
          }
          cur_shared_idx[k] = 0;
          offA -= (plan.sdim[k] - 1) * shared_stride_a[k];
          offB -= (plan.sdim[k] - 1) * shared_stride_b[k];
        }
      }
      out.data[r] = acc;
      for (int t = plan.RR - 1; t >= 0; --t) { // step the result odometer
        if (++cur_result_idx[t] < plan.rdim[t]) {
          baseA += result_stride_a[t];
          baseB += result_stride_b[t];
          break;
        }
        cur_result_idx[t] = 0;
        baseA -= (plan.rdim[t] - 1) * result_stride_a[t];
        baseB -= (plan.rdim[t] - 1) * result_stride_b[t];
      }
    }
  }

  /// @brief Copy a view into an owned @ref DynTensor (used for the single-operand fold base case
  ///        and for @ref scale — both produce a result that must outlive its sibling operands).
  constexpr DynTensor own(detail::View v)
  {
    DynTensor out;
    out.rank = v.rank;
    std::size_t n = 1;
    for (int i = 0; i < v.rank; ++i) {
      out.ids[i] = v.ids[i];
      out.dims[i] = v.dims[i];
      n *= v.dims[i];
    }
    out.data.assign(v.data, v.data + n);
    return out;
  }

  /// @brief Contract a single argument — base case of the variadic fold.
  template <class T> constexpr DynTensor contract_all(const T &x) { return own(view_of(x)); }
  /// @brief Left-fold a contraction chain (closing a loop yields a scalar trace). Each argument
  ///        may be a typed builder @ref DTensor or a @ref DynTensor from a nested fold; operands
  ///        are viewed in place and the running intermediate ping-pongs between two reusable
  ///        buffers, so heap traffic is O(1) and the big intermediates never sit on the stack.
  template <class A, class B, class... Rest>
  constexpr DynTensor contract_all(const A &a, const B &b, const Rest &...rest)
  {
    const detail::View vs[] = {view_of(a), view_of(b), view_of(rest)...};
    constexpr std::size_t N = 2 + sizeof...(Rest);
    DynTensor ping, pong;
    DynTensor *cur = &ping, *oth = &pong;
    detail::View acc = vs[0];
    for (std::size_t i = 1; i < N; ++i) {
      contract_into(*cur, acc, vs[i]); // reads acc (in *oth or vs[0]), writes *cur (distinct)
      acc = view_of(*cur);
      DynTensor *tmp = cur;
      cur = oth;
      oth = tmp;
    }
    return std::move(*oth); // last write landed in *oth after the final swap
  }

  /// @brief Scale a typed builder or a nested fold result by a compile-time complex factor.
  template <Cx C, class T> constexpr DynTensor scale(const T &x)
  {
    DynTensor out = own(view_of(x));
    for (auto &z : out.data)
      z = C * z;
    return out;
  }

  // ---- entrywise sum (EAddPlan reindexing; numeric add) ------------------------

  /// @brief Sum views `a`,`b` (B reindexed into A's axis order by identity, via
  ///        @ref numtracer::core::EAddPlan), writing into `out` (capacity reused; `out` must not
  ///        alias either operand) — the eager structure-summation primitive.
  constexpr void add_into(DynTensor &out, detail::View a, detail::View b)
  {
    std::array<int, core::kMaxAxisRank> a_ids{}, a_dims{}, b_ids{}, b_dims{};
    for (int i = 0; i < a.rank; ++i) {
      a_ids[i] = a.ids[i];
      a_dims[i] = a.dims[i];
    }
    for (int i = 0; i < b.rank; ++i) {
      b_ids[i] = b.ids[i];
      b_dims[i] = b.dims[i];
    }
    const core::EAddPlan plan = core::make_eaddplan(a_ids, a_dims, a.rank, b_ids, b_dims);
    out.rank = a.rank;
    for (int i = 0; i < a.rank; ++i) {
      out.ids[i] = a.ids[i];
      out.dims[i] = a.dims[i];
    }
    out.data.resize(plan.n);
    for (std::size_t r = 0; r < plan.n; ++r)
      out.data[r] = a.data[r] + b.data[core::add_bindex(r, plan)];
  }

  /// @brief Sum a single argument — base case of the variadic fold.
  template <class T> constexpr DynTensor add_all(const T &x) { return own(view_of(x)); }
  /// @brief Left-fold an entrywise sum over same-axis tensors (the eager vertex-structure sum).
  ///        Same ping-pong/view discipline as @ref contract_all.
  template <class A, class B, class... Rest> constexpr DynTensor add_all(const A &a, const B &b, const Rest &...rest)
  {
    const detail::View vs[] = {view_of(a), view_of(b), view_of(rest)...};
    constexpr std::size_t N = 2 + sizeof...(Rest);
    DynTensor ping, pong;
    DynTensor *cur = &ping, *oth = &pong;
    detail::View acc = vs[0];
    for (std::size_t i = 1; i < N; ++i) {
      add_into(*cur, acc, vs[i]);
      acc = view_of(*cur);
      DynTensor *tmp = cur;
      cur = oth;
      oth = tmp;
    }
    return std::move(*oth);
  }

NT_END_NO_LOOP_VECTORIZE

  // ============================================================================
  // Sector builders. Constant builders are constexpr/nullary (colour, gamma) so a renv-free
  // component folds to a compile-time constant; momentum builders take the same `double renv[]`
  // the codegen fills.
  // ============================================================================

  // ---- Lorentz ---------------------------------------------------------------

  /// @brief The Euclidean metric @f$\delta_{\mu\nu}@f$.
  template <int Mu, int Nu> constexpr auto metric()
  {
    DTensor<AxList<Ax<Mu, 4>, Ax<Nu, 4>>> t{};
    for (int i = 0; i < 4; ++i)
      t.data[i * 4 + i] = Cx{1, 0};
    return t;
  }

  /// @brief A frame-fixed momentum 4-vector. `Mask` flags the structurally-nonzero Lorentz
  ///        components; `renv` stores only those, packed: the k-th present component (in
  ///        ascending Lorentz order) is `renv[Base+k]`.
  template <int Lbl, int Base, int Mask> auto vector(const double *renv)
  {
    DTensor<AxList<Ax<Lbl, 4>>> t{};
    int k = 0;
    for (int i = 0; i < 4; ++i)
      t.data[i] = ((Mask >> i) & 1) ? Cx{renv[Base + k++], 0} : Cx{0, 0};
    return t;
  }

  /// @brief The transverse projector @f$P_{\mu\nu}(l) = \delta_{\mu\nu} - l_\mu l_\nu/l^2@f$
  ///        (`renv[Inv]` holds @f$1/l^2@f$; masked components are zero, present ones packed
  ///        from `renv[Base]` in Lorentz order).
  template <int Mu, int Nu, int Base, int Mask, int Inv> auto transverse_projector(const double *renv)
  {
    DTensor<AxList<Ax<Mu, 4>, Ax<Nu, 4>>> t{};
    double q[4];
    int k = 0;
    for (int i = 0; i < 4; ++i)
      q[i] = ((Mask >> i) & 1) ? renv[Base + k++] : 0.0;
    const double inv = renv[Inv];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        t.data[i * 4 + j] = Cx{(i == j ? 1.0 : 0.0) - q[i] * q[j] * inv, 0};
    return t;
  }

  /// @brief The longitudinal projector @f$P^L_{\mu\nu}(l) = l_\mu l_\nu/l^2@f$.
  template <int Mu, int Nu, int Base, int Mask, int Inv> auto longitudinal_projector(const double *renv)
  {
    DTensor<AxList<Ax<Mu, 4>, Ax<Nu, 4>>> t{};
    double q[4];
    int k = 0;
    for (int i = 0; i < 4; ++i)
      q[i] = ((Mask >> i) & 1) ? renv[Base + k++] : 0.0;
    const double inv = renv[Inv];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        t.data[i * 4 + j] = Cx{q[i] * q[j] * inv, 0};
    return t;
  }

  // ---- Dirac (constant) ------------------------------------------------------

  /// @brief @f$\gamma^\mu@f$ carrying a free Lorentz axis (rank-3: Lorentz ⊗ two spinor).
  template <int MuLbl, int Din, int Dout> constexpr auto gamma_axis()
  {
    DTensor<AxList<Ax<MuLbl, 4>, Ax<Din, 4>, Ax<Dout, 4>>> t{};
    for (int mu = 0; mu < 4; ++mu)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          t.data[(mu * 4 + i) * 4 + j] = dirac::gamma_entry(mu, i, j);
    return t;
  }

  /// @brief The spinor identity (contract with `M` to take @f$\mathrm{tr}(M)@f$).
  template <int Din, int Dout> constexpr auto identity()
  {
    DTensor<AxList<Ax<Din, 4>, Ax<Dout, 4>>> t{};
    for (int i = 0; i < 4; ++i)
      t.data[i * 4 + i] = Cx{1, 0};
    return t;
  }

  // ---- SU(N) colour (constant; read straight from the typed-out tables) ------

  /// @brief The generator @f$(T^A)_{IJ}@f$ on one adjoint and two fundamental axes.
  template <int N, int A, int I, int J> constexpr auto T()
  {
    constexpr int D = N * N - 1;
    DTensor<AxList<Ax<A, D>, Ax<I, N>, Ax<J, N>>> t{};
    for (int a = 0; a < D; ++a)
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
          const auto z = sun::SUNData<N>::generators[a].data[i * N + j];
          t.data[(a * N + i) * N + j] = Cx{z.real(), z.imag()};
        }
    return t;
  }

  /// @brief The structure constant @f$f^{ABC}@f$ on three adjoint axes.
  template <int N, int A, int B, int C> constexpr auto f()
  {
    constexpr int D = N * N - 1;
    DTensor<AxList<Ax<A, D>, Ax<B, D>, Ax<C, D>>> t{};
    for (const auto &e : sun::SUNData<N>::f_nonzeros)
      t.data[(e.a * D + e.b) * D + e.c] = Cx{e.v, 0};
    return t;
  }

  /// @brief The Kronecker delta on two fundamental axes.
  template <int N, int I, int J> constexpr auto delta_fund()
  {
    DTensor<AxList<Ax<I, N>, Ax<J, N>>> t{};
    for (int i = 0; i < N; ++i)
      t.data[i * N + i] = Cx{1, 0};
    return t;
  }

  /// @brief The Kronecker delta on two adjoint axes.
  template <int N, int A, int B> constexpr auto delta_adj()
  {
    constexpr int D = N * N - 1;
    DTensor<AxList<Ax<A, D>, Ax<B, D>>> t{};
    for (int a = 0; a < D; ++a)
      t.data[a * D + a] = Cx{1, 0};
    return t;
  }

} // namespace numtracer::dense
