/// @file cmat.hpp
/// @brief A small compile-time-sized dense matrix, generic in its entry type.
///
/// This header is a foundational leaf — it has no dependency on the Lorentz/Dirac/SU(N) sector
/// headers. `Matrix<N>` (complex entries) holds the typed-out SU(N) generators in `sun/sun_data.hpp`;
/// `Matrix<4, Poly>` is the spinor matrix of the numeric Dirac trace (`engine/spinor_mat.hpp`).
/// The type is an aggregate, so it stays brace-initialisable and usable in `constexpr` tables; the
/// product and trace below are generic in both `N` and the entry type.
#pragma once

#include <array>
#include <complex>
#include <utility> // std::move

namespace numtracer
{

  /// @brief A compile-time-sized N×N matrix stored on the stack (row-major).
  /// @tparam N The matrix dimension.
  /// @tparam T The entry type (complex by default).
  template <int N, typename T = std::complex<double>> struct Matrix {
    std::array<T, N * N> data{}; ///< Row-major element storage.
    /// @brief Mutable element access.
    /// @param i Row index.
    /// @param j Column index.
    /// @return Reference to element `(i,j)`.
    constexpr T &operator()(int i, int j) { return data[i * N + j]; }
    /// @brief Const element access.
    /// @param i Row index.
    /// @param j Column index.
    /// @return Reference to element `(i,j)`.
    constexpr const T &operator()(int i, int j) const { return data[i * N + j]; }
  };

  /// @brief Matrix product `A · B`.
  ///
  /// Each entry starts from its first product rather than from a zero, so `T` needs only `+` and `*`
  /// (a `Poly` zero would need the frame's symbol count).
  template <int N, typename T> Matrix<N, T> matmul(const Matrix<N, T> &A, const Matrix<N, T> &B)
  {
    static_assert(N >= 1, "matmul needs N >= 1");
    Matrix<N, T> C;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) {
        T s = A(i, 0) * B(0, j);
        for (int k = 1; k < N; ++k)
          s = s + A(i, k) * B(k, j);
        C(i, j) = std::move(s);
      }
    return C;
  }

  /// @brief Matrix trace `Σ_i A(i,i)`.
  template <int N, typename T> T mtrace(const Matrix<N, T> &A)
  {
    static_assert(N >= 1, "mtrace needs N >= 1");
    T s = A(0, 0);
    for (int i = 1; i < N; ++i)
      s = s + A(i, i);
    return s;
  }

} // namespace numtracer
