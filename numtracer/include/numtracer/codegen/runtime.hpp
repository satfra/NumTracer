/// @file runtime.hpp
/// @brief Minimal runtime-support symbols for NumTracer-generated kernels.
///
/// A generated kernel evaluates straight-line real/complex arithmetic and a handful of
/// scalar math helpers (`powr`, `pow`, `sqrt`, `fma`) plus a runtime `complex` type. The
/// codegen leaves these unqualified and brings them into scope with
/// `using namespace <SupportNamespace>;` / `using namespace <SupportNamespace>::compute;`.
///
/// This header is the **default** support namespace (`numtracer`), so default-generated
/// kernels are self-contained against the NumTracer headers alone. A consumer that already
/// provides an equivalent API (its own `complex` alias and `compute::powr/...`) can point the
/// codegen at it instead via the `"RuntimeInclude"` / `"SupportNamespace"` options of
/// `MakeNTKernel` — the generated code itself is agnostic to which one supplies the names.
#pragma once

#include <cmath>
#include <complex>

namespace numtracer
{
  /// Runtime complex scalar used by generated kernels (distinct from the compile-time
  /// non-type-parameter complex `numtracer::Cx`). Left as a template alias so a device/CUDA
  /// support header can substitute a device-safe complex without touching the emitted code.
  template <class T> using complex = std::complex<T>;

  namespace compute
  {
    /// @brief Compile-time integer power, x^n (handles negative/zero exponents).
    /// @tparam n Integer exponent.
    template <int n, typename NumberType>
      requires requires(NumberType x) {
        x * x;
        NumberType(1.) / x;
      }
    constexpr inline NumberType powr(const NumberType x)
    {
      if constexpr (n == 0)
        return NumberType(1.);
      else if constexpr (n < 0)
        return NumberType(1.) / powr<-n, NumberType>(x);
      else if constexpr (n == 1)
        return x;
      else if constexpr (n % 2 == 0)
        return powr<n / 2>(x) * powr<n / 2>(x);
      else
        return powr<n / 2>(x) * powr<n / 2>(x) * x;
    }

    using std::fma;
    using std::pow;
    using std::sqrt;

    /// Complex fused multiply-add: a kernel whose projector carries an `i` assembles a complex
    /// intermediate via `fma()`, and there is no `std::fma` for `std::complex`. The real
    /// `fma(double,double,double)` still binds `std::fma` exactly, so real kernels are unaffected.
    inline std::complex<double> fma(const std::complex<double> &a, const std::complex<double> &b,
                                    const std::complex<double> &c)
    {
      return a * b + c;
    }
  } // namespace compute
} // namespace numtracer
