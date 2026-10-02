/// @file sun_data.hpp
/// @brief Compile-time SU(2)/SU(3) colour tables (generators and structure constants).
///
/// The generators @f$T^a@f$ and structure constants @f$f^{abc}@f$ depend only on
/// `N`, so for the physically relevant cases they are not built at runtime: they
/// live here as `static constexpr` numerical literals. `std::complex`'s
/// constructor is `constexpr` in C++20, so the generator tables are genuine
/// compile-time constants, and the @f$f^{abc}@f$ are plain reals — no `constexpr`
/// complex arithmetic or `sqrt` is needed.
///
/// These values were generated once from the generalized-Gell-Mann construction at full
/// (`%.17g`) precision; `tests/test_sun_tables.cpp` cross-checks them against that construction
/// (`build_oracle` in `network/sun_net.hpp`, which also serves `N >= 4`), so they cannot silently
/// drift.
///
/// Conventions: SU(2) uses @f$T^a = \sigma^a/2@f$ (so @f$f^{abc} =
/// \epsilon^{abc}@f$); SU(3) uses the Gell-Mann @f$\lambda^a/2@f$, with the
/// @f$\lambda_8@f$ diagonal at @f$\pm 1/(2\sqrt{3})@f$ and the @f$\sqrt{3}/2@f$
/// entries of @f$f@f$ equal to `0.8660254037844386`.
#pragma once

#include "numtracer/core/cmat.hpp" // Mat<N>, matmul, trace (the dense complex-matrix leaf)
#include "numtracer/core/cx.hpp"   // Cx (constexpr complex), approx
#include "numtracer/core/config.hpp" // NT_THROW (exception-optional guard for -fno-exceptions builds)

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace numtracer::sun {

/// @brief Internal helpers shared by the SU(N) builder and its typed-out tables.
namespace sun_detail {

/// @brief One nonzero adjoint structure constant @f$f^{abc} = v@f$.
///
/// Used by both the runtime builder (`network/sun_net.hpp`) and the typed-out `constexpr` tables
/// below.
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

} // namespace numtracer::sun
