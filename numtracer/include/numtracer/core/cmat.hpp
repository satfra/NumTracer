/// @file cmat.hpp
/// @brief A small compile-time-sized dense complex matrix and its basic operations.
///
/// This header is a foundational leaf — it has no dependency on the Lorentz/Dirac/SU(N)
/// sector headers. It provides a stack-allocated `N×N` complex matrix (`std::complex<double>`
/// entries) with element access, a sparse-aware matrix product, and a trace. The SU(N) sector
/// (`sun/sun_data.hpp`, `network/sun_net.hpp`) builds its generators and structure constants on
/// top of these primitives; the matmul skips structurally-zero left entries because those
/// generators are sparse.
#pragma once

#include <array>
#include <complex>

namespace numtracer {

/// @brief A compile-time-sized N×N complex matrix stored on the stack (row-major).
/// @tparam N The matrix dimension.
template <int N> struct Mat {
  std::array<std::complex<double>, N * N> data{}; ///< Row-major element storage.
  /// @brief Mutable element access.
  /// @param i Row index.
  /// @param j Column index.
  /// @return Reference to element `(i,j)`.
  std::complex<double> &operator()(int i, int j) { return data[i * N + j]; }
  /// @brief Const element access.
  /// @param i Row index.
  /// @param j Column index.
  /// @return The value of element `(i,j)`.
  std::complex<double> operator()(int i, int j) const { return data[i * N + j]; }
};

/// @brief Matrix product, skipping structurally zero entries of the left operand.
///
/// The generators are sparse, so skipping zero `x(i,k)` saves work.
/// @tparam N The matrix dimension.
/// @param x The left matrix.
/// @param y The right matrix.
/// @return The product `x * y`.
template <int N> Mat<N> matmul(const Mat<N> &x, const Mat<N> &y) {
  Mat<N> r;
  for (int i = 0; i < N; ++i)
    for (int k = 0; k < N; ++k) {
      const std::complex<double> xik = x(i, k);
      if (xik == std::complex<double>{0, 0})
        continue; // generators are sparse — skip structural zeros
      for (int j = 0; j < N; ++j)
        r(i, j) += xik * y(k, j);
    }
  return r;
}

/// @brief Matrix trace.
/// @tparam N The matrix dimension.
/// @param x The matrix.
/// @return @f$\sum_i x_{ii}@f$.
template <int N> std::complex<double> trace(const Mat<N> &x) {
  std::complex<double> t{0, 0};
  for (int i = 0; i < N; ++i)
    t += x(i, i);
  return t;
}

} // namespace numtracer
