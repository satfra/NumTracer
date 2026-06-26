/// @file sun_data.hpp
/// @brief Compile-time SU(2)/SU(3) colour tables plus the @ref numtracer::sun::SUNBuilder
///        sector class (runtime oracle for the generators + structure constants).
///
/// The generators @f$T^a@f$ and structure constants @f$f^{abc}@f$ depend only on
/// `N`, so for the physically relevant cases they are not built at runtime: they
/// live here as `static constexpr` numerical literals. `std::complex`'s
/// constructor is `constexpr` in C++20, so the generator tables are genuine
/// compile-time constants, and the @f$f^{abc}@f$ are plain reals — no `constexpr`
/// complex arithmetic or `sqrt` is needed.
///
/// These values were generated once from the runtime @ref numtracer::sun::SUNBuilder
/// at full (`%.17g`) precision; `tests/test_core.cpp` cross-checks them against
/// that builder (the oracle), so they can never silently drift. The runtime
/// builder is retained for `N >= 4` and as that oracle.
///
/// Conventions: SU(2) uses @f$T^a = \sigma^a/2@f$ (so @f$f^{abc} =
/// \epsilon^{abc}@f$); SU(3) uses the Gell-Mann @f$\lambda^a/2@f$, with the
/// @f$\lambda_8@f$ diagonal at @f$\pm 1/(2\sqrt{3})@f$ and the @f$\sqrt{3}/2@f$
/// entries of @f$f@f$ equal to `0.8660254037844386`.
#pragma once

#include "numtracer/core/cmat.hpp" // Mat<N>, matmul, trace (the dense complex-matrix leaf)
#include "numtracer/core/cx.hpp"   // Cx (constexpr complex), approx

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace numtracer::sun {

/// @brief Internal helpers shared by the SU(N) oracle and its typed-out tables.
namespace sun_detail {

/// @brief One nonzero adjoint structure constant @f$f^{abc} = v@f$.
///
/// Used by both the runtime @ref numtracer::sun::SUNBuilder store and the typed-out
/// `constexpr` tables below.
struct FEntry {
  int a;    ///< First adjoint index.
  int b;    ///< Second adjoint index.
  int c;    ///< Third adjoint index.
  double v; ///< The structure-constant value @f$f^{abc}@f$.
};

} // namespace sun_detail

/// @brief Typed-out compile-time SU(N) data; only specialized for tabulated `N`.
/// @tparam N The colour group rank.
template <int N> struct SUNData;          // only specialized for the typed-out N
/// @brief Whether a @ref numtracer::sun::SUNData specialization exists for `N`.
/// @tparam N The colour group rank.
template <int N> inline constexpr bool kHasSUNData = false;

/// @brief Compile-time SU(2) data: the 3 generators and the nonzero @f$f^{abc}@f$.
template <> struct SUNData<2> {
  /// @brief The 3 fundamental generators @f$T^a = \sigma^a/2@f$.
  static constexpr std::array<Mat<2>, 3> generators = {{
    Mat<2>{{{std::complex<double>{0,0}, std::complex<double>{0.5,0}, std::complex<double>{0.5,0}, std::complex<double>{0,0}}}},
    Mat<2>{{{std::complex<double>{0,0}, std::complex<double>{0,-0.5}, std::complex<double>{0,0.5}, std::complex<double>{0,0}}}},
    Mat<2>{{{std::complex<double>{0.5,0}, std::complex<double>{0,0}, std::complex<double>{0,0}, std::complex<double>{-0.5,0}}}},
  }};
  /// @brief The nonzero structure constants @f$f^{abc} = \epsilon^{abc}@f$.
  static constexpr std::array<sun_detail::FEntry, 6> f_nonzeros = {{
    {0,1,2,1}, {0,2,1,-1},
    {1,0,2,-1}, {1,2,0,1},
    {2,0,1,1}, {2,1,0,-1},
  }};
};
template <> inline constexpr bool kHasSUNData<2> = true; ///< SU(2) data is tabulated.

/// @brief Compile-time SU(3) data: the 8 generators and the nonzero @f$f^{abc}@f$.
template <> struct SUNData<3> {
  /// @brief The 8 fundamental generators @f$T^a = \lambda^a/2@f$ (Gell-Mann).
  static constexpr std::array<Mat<3>, 8> generators = {{
    Mat<3>{{{std::complex<double>{0,0},std::complex<double>{0.5,0},std::complex<double>{0,0},std::complex<double>{0.5,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0}}}},
    Mat<3>{{{std::complex<double>{0,0},std::complex<double>{0,-0.5},std::complex<double>{0,0},std::complex<double>{0,0.5},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0}}}},
    Mat<3>{{{std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0.5,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0.5,0},std::complex<double>{0,0},std::complex<double>{0,0}}}},
    Mat<3>{{{std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,-0.5},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0.5},std::complex<double>{0,0},std::complex<double>{0,0}}}},
    Mat<3>{{{std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0.5,0},std::complex<double>{0,0},std::complex<double>{0.5,0},std::complex<double>{0,0}}}},
    Mat<3>{{{std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,-0.5},std::complex<double>{0,0},std::complex<double>{0,0.5},std::complex<double>{0,0}}}},
    Mat<3>{{{std::complex<double>{0.5,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{-0.5,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0}}}},
    Mat<3>{{{std::complex<double>{0.28867513459481287,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0.28867513459481287,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{0,0},std::complex<double>{-0.57735026918962573,0}}}},
  }};
  /// @brief The 54 nonzero structure constants @f$f^{abc}@f$.
  static constexpr std::array<sun_detail::FEntry, 54> f_nonzeros = {{
    {0,1,6,1}, {0,2,5,0.5}, {0,3,4,-0.5}, {0,4,3,0.5}, {0,5,2,-0.5}, {0,6,1,-1},
    {1,0,6,-1}, {1,2,4,0.5}, {1,3,5,0.5}, {1,4,2,-0.5}, {1,5,3,-0.5}, {1,6,0,1},
    {2,0,5,-0.5}, {2,1,4,-0.5}, {2,3,6,0.5}, {2,3,7,0.8660254037844386}, {2,4,1,0.5},
    {2,5,0,0.5}, {2,6,3,-0.5}, {2,7,3,-0.8660254037844386},
    {3,0,4,0.5}, {3,1,5,-0.5}, {3,2,6,-0.5}, {3,2,7,-0.8660254037844386}, {3,4,0,-0.5},
    {3,5,1,0.5}, {3,6,2,0.5}, {3,7,2,0.8660254037844386},
    {4,0,3,-0.5}, {4,1,2,0.5}, {4,2,1,-0.5}, {4,3,0,0.5}, {4,5,6,-0.5},
    {4,5,7,0.8660254037844386}, {4,6,5,0.5}, {4,7,5,-0.8660254037844386},
    {5,0,2,0.5}, {5,1,3,0.5}, {5,2,0,-0.5}, {5,3,1,-0.5}, {5,4,6,0.5},
    {5,4,7,-0.8660254037844386}, {5,6,4,-0.5}, {5,7,4,0.8660254037844386},
    {6,0,1,1}, {6,1,0,-1}, {6,2,3,0.5}, {6,3,2,-0.5}, {6,4,5,-0.5}, {6,5,4,0.5},
    {7,2,3,0.8660254037844386}, {7,3,2,-0.8660254037844386}, {7,4,5,0.8660254037844386}, {7,5,4,-0.8660254037844386},
  }};
};
template <> inline constexpr bool kHasSUNData<3> = true; ///< SU(3) data is tabulated.

namespace sun_detail {
/// @brief Compile-time lookup of a structure constant in the @ref numtracer::sun::SUNData table.
/// @tparam N The colour group rank.
/// @tparam a First adjoint index.
/// @tparam b Second adjoint index.
/// @tparam c Third adjoint index.
/// @return @f$f^{abc}@f$ as a @ref numtracer::Cx (zero if absent from the table).
template <int N, int a, int b, int c> constexpr Cx f_val() {
  for (const auto &e : SUNData<N>::f_nonzeros)
    if (e.a == a && e.b == b && e.c == c) return Cx{e.v, 0.0};
  return Cx{0.0, 0.0};
}
/// @brief Compile-time lookup of a generator entry @f$(T^a)_{ij}@f$.
/// @tparam N The colour group rank.
/// @tparam a The adjoint (generator) index.
/// @tparam i The row (fundamental) index.
/// @tparam j The column (fundamental) index.
/// @return @f$(T^a)_{ij}@f$ as a @ref numtracer::Cx.
template <int N, int a, int i, int j> constexpr Cx gen_val() {
  const auto z = SUNData<N>::generators[a].data[i * N + j];
  return Cx{z.real(), z.imag()};
}
} // namespace sun_detail

/// @brief The SU(N) sector: a runtime numeric oracle for the colour algebra.
///
/// The constructor builds the fundamental generators (generalized Gell-Mann) and the
/// structure constants @f$f^{abc} = -2 i\,\mathrm{tr}([T^a,T^b]T^c)@f$. This is the
/// source of truth the typed-out @ref SUNData tables are validated against, and the
/// substrate the numeric colour contraction (`network/sun_net.hpp`) mirrors.
/// @tparam N The colour group rank.
template <int N> class SUNBuilder {
public:
  /// @brief The adjoint dimension @f$N^2-1@f$ (number of generators).
  static constexpr int kAdjDim = N * N - 1;

  /// @brief Construct the oracle, building the generators and structure constants.
  SUNBuilder() {
    build_generators();
    build_structure_constants();
  }
  /// @brief The number of colours `N`.
  /// @return `N`.
  static constexpr int num_colors() { return N; }
  /// @brief The adjoint dimension.
  /// @return @f$N^2-1@f$.
  static constexpr int adj_dim() { return kAdjDim; }
  /// @brief Access a fundamental generator built by the oracle.
  /// @param a The adjoint (generator) index.
  /// @return The matrix @f$T^a@f$.
  const Mat<N> &gen(int a) const { return gens_[a]; }                 // fundamental generator T^a
  /// @brief The nonzero structure constants built by the oracle.
  /// @return The list of @ref numtracer::sun_detail::FEntry records.
  const std::vector<sun_detail::FEntry> &f_nonzeros() const { return f_nz_; }

private:
  std::array<Mat<N>, kAdjDim> gens_{}; ///< The oracle-built fundamental generators.
  std::vector<sun_detail::FEntry> f_nz_; ///< The oracle-built nonzero structure constants.

  /// @brief Build the fundamental generators (generalized Gell-Mann) into @ref gens_.
  /// @throws std::logic_error if the generator count does not match @ref kAdjDim.
  void build_generators() {
    const std::complex<double> I{0, 1};
    int idx = 0;
    for (int k = 0; k < N; ++k)
      for (int l = k + 1; l < N; ++l) {
        Mat<N> m;
        m(k, l) = std::complex<double>{0.5, 0};
        m(l, k) = std::complex<double>{0.5, 0};
        gens_[idx++] = m;
        Mat<N> n;
        n(k, l) = -I * 0.5;
        n(l, k) = I * 0.5;
        gens_[idx++] = n;
      }
    for (int l = 1; l <= N - 1; ++l) {
      Mat<N> m;
      const double norm = std::sqrt(1.0 / (2.0 * l * (l + 1)));
      for (int j = 0; j < l; ++j) m(j, j) = std::complex<double>{norm, 0};
      m(l, l) = std::complex<double>{-norm * l, 0};
      gens_[idx++] = m;
    }
    if (idx != kAdjDim) throw std::logic_error("SUNBuilder: wrong generator count");
  }
  /// @brief Build the nonzero structure constants @f$f^{abc} = -2 i\,\mathrm{tr}([T^a,T^b]T^c)@f$.
  ///
  /// Stores every entry whose magnitude exceeds a `1e-12` round-off tolerance into @ref f_nz_.
  void build_structure_constants() {
    // entries below this magnitude are floating-point round-off in the commutator trace, not genuine
    // nonzero structure constants, so they are dropped.
    constexpr double kStructConstTol = 1e-12;
    for (int a = 0; a < kAdjDim; ++a)
      for (int b = 0; b < kAdjDim; ++b) {
        const Mat<N> ab = matmul(gens_[a], gens_[b]);
        const Mat<N> ba = matmul(gens_[b], gens_[a]);
        Mat<N> comm;
        for (std::size_t i = 0; i < comm.data.size(); ++i) comm.data[i] = ab.data[i] - ba.data[i];
        for (int c = 0; c < kAdjDim; ++c) {
          const std::complex<double> tr = trace(matmul(comm, gens_[c]));
          const double val = (std::complex<double>{0, -2} * tr).real();
          if (std::fabs(val) > kStructConstTol) f_nz_.push_back({a, b, c, val});
        }
      }
  }
};

} // namespace numtracer::sun
