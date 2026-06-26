/// @file cx.hpp
/// @brief Minimal `constexpr` complex number usable as a non-type template
///        parameter.
///
/// `std::complex` arithmetic is not fully `constexpr` before C++23, so the
/// expression algebra and the tensor entries use this minimal complex type
/// instead. It is a *structural type* (all members are public and of literal
/// type), which C++20 lets us pass as a non-type template parameter — that is
/// what makes a compile-time complex literal such as `expr::Lit<Cx{0,1}>`
/// possible. This header is foundational and intentionally has no dependency on
/// the Lorentz/Dirac/SU(N) sector headers.
#pragma once

namespace numtracer
{

  /// @brief A complex number with `double` components, fully usable at compile time.
  ///
  /// Defined as a structural type so a `Cx` value can appear as a non-type template
  /// parameter (see @ref numtracer::expr::Lit). Only addition and multiplication are
  /// provided, which is all the expression algebra needs.
  struct Cx {
    double re = 0; ///< Real part.
    double im = 0; ///< Imaginary part.
    /// @brief Complex addition.
    /// @param o The right-hand operand.
    /// @return The sum `(re+o.re, im+o.im)`.
    constexpr Cx operator+(Cx o) const { return {re + o.re, im + o.im}; }
    /// @brief Complex multiplication.
    /// @param o The right-hand operand.
    /// @return The product `(re*o.re - im*o.im, re*o.im + im*o.re)`.
    constexpr Cx operator*(Cx o) const { return {re * o.re - im * o.im, re * o.im + im * o.re}; }
  };

  /// @brief Absolute value of a real number (a `constexpr` `std::abs` for `double`).
  /// @param x The input value.
  /// @return `x` if `x >= 0`, otherwise `-x`.
  constexpr double cabs_real(double x) { return x < 0 ? -x : x; }

  /// @brief Test whether two complex numbers agree to within a tolerance.
  /// @param a First value.
  /// @param b Second value.
  /// @param tol Maximum allowed absolute difference per component (default `1e-9`).
  /// @return `true` if both the real and imaginary parts differ by at most `tol`.
  constexpr bool approx(Cx a, Cx b, double tol = 1e-9)
  {
    return cabs_real(a.re - b.re) <= tol && cabs_real(a.im - b.im) <= tol;
  }

} // namespace numtracer
