// NumTracer — standalone shim so the copied FormTracer kernels compile and run
// with NO DiFfRG dependency. It provides:
//   * the DiFfRG / DiFfRG::compute names the generated kernels expect
//     (powr, pow, sqrt, fma, the interpolator/coordinate/memory type stubs,
//      KOKKOS_FORCEINLINE_FUNCTION),
//   * a fixed analytic Regulator and a fixed set of dressing functions used
//     IDENTICALLY by both the reference kernel and the NumTracer assembly.
//
// The specific regulator/dressing formulas are arbitrary but fixed: the test
// checks that the two evaluation PATHS agree as algebraic functions, so any
// consistent definitions used on both sides suffice.
#pragma once

#include <array>
#include <cmath>
#include <complex>

#define KOKKOS_FORCEINLINE_FUNCTION inline

namespace DiFfRG
{

  // The complex type the generated INV kernels resolve unqualified (they are emitted without
  // `using std::complex;` so that a CUDA shim can substitute cuda::std::complex — std::complex
  // arithmetic silently miscompiles to 0 in nvcc device code). Here it is plain std::complex.
  template <class T> using complex = std::complex<T>;

  // ---- math names used by generated kernels ----------------------------------
  namespace compute
  {
    /**
     * @brief A compile-time evaluatable power function for whole number exponents
     *
     * @tparam n Exponent of type int
     * @tparam RF Type of argument
     * @param x Argument
     * @return x^n
     */
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
    // complex fused multiply-add: the dense kernel of a Dirac VERTEX (whose projector carries an
    // `i`) assembles a complex integrand via fma(); there is no std::fma for std::complex. Real
    // fma(double,double,double) still binds std::fma exactly, so real kernels are unaffected.
    inline std::complex<double> fma(const std::complex<double> &a, const std::complex<double> &b,
                                    const std::complex<double> &c)
    {
      return a * b + c;
    }
  } // namespace compute

  // ---- type stubs the generated kernel signatures mention --------------------
  struct GPU_memory {
  };
  template <class T> struct LogarithmicCoordinates1D {
  };

  // A 1-D "interpolator" that just wraps a function pointer (double -> double).
  // Aliased to the generated SplineInterpolator1D type name below.
  struct Fn {
    double (*f)(double) = nullptr;
    double operator()(double x) const { return f(x); }
  };
  template <class A, class B, class C> using SplineInterpolator1D = Fn;

  // Per-component (group-diagonal) dressing access: a kernel passes a std::array of interpolators
  // (one per colour/flavour component) and the generated integrand calls ntDiagDress(arr, a, scale)
  // for the surviving component `a` of the SU(N) fold. Indexed lookup + evaluation at the scale.
  template <class Arr> inline auto ntDiagDress(const Arr &arr, int a, double scale)
  {
    return arr[static_cast<std::size_t>(a)](scale);
  }

  // ---- a fixed analytic regulator (same on both paths) -----------------------
  // Exponential-type: RB(k2,p2) = k2 * exp(-p2/k2). RBdot is a fixed companion.
  // The fermionic RB/RFdot are nonzero analytic companions (the exact form is
  // arbitrary but FIXED and identical on both paths — the Zq quark loop needs them
  // live so the regularized quark propagator and its regulator-dot are exercised).
  struct ShimRegulator {
    static double RB(double k2, double p2) { return k2 * std::exp(-p2 / k2); }
    static double RBdot(double k2, double p2) { return 2.0 * (k2 + p2) * std::exp(-p2 / k2); }
    static double RF(double k2, double p2) { return std::sqrt(k2) * std::exp(-p2 / k2); }
    static double RFdot(double k2, double p2) { return (k2 + 2.0 * p2) / std::sqrt(k2) * std::exp(-p2 / k2); }
    static double dq2RB(double, double) { return 0.0; }
    static double dq2RF(double, double) { return 0.0; }
  };

} // namespace DiFfRG

// ---- fixed dressing functions, shared by both paths ------------------------
// Distinct smooth positive functions of one argument.
namespace shimdress
{
  inline double ZA(double x) { return 1.0 / (1.0 + 0.5 * x * x) + 0.7; }
  inline double Zc(double x) { return 1.0 / (1.0 + 0.3 * x) + 0.9; }
  inline double ZAcbc(double x) { return 1.2 + 0.4 * std::sin(0.5 * x); }
  inline double ZA3(double x) { return 1.1 + 0.2 * std::cos(0.3 * x); }
  inline double ZA4(double x) { return 1.05 + 0.1 * std::sin(0.7 * x); }
  inline double dtZA(double x) { return -0.2 / (1.0 + x * x); }
  inline double dtZc(double x) { return -0.15 / (1.0 + 0.2 * x); }
  // quark sector: quark-gluon vertex, kinetic dressing Zq, mass Mq (kept nonzero so
  // the +Mq^2 term in the quark propagator denominator is genuinely exercised), dtZq.
  inline double ZAqbq1(double x) { return 1.15 + 0.25 * std::cos(0.4 * x); }
  // quark-gluon vertex tensor structures 4 (p^mu) and 7 (sigma^{mu nu}) of the {1,4,7} basis.
  inline double ZAqbq4(double x) { return 0.6 + 0.2 * std::sin(0.35 * x); }
  inline double ZAqbq7(double x) { return 0.4 + 0.15 * std::cos(0.55 * x); }
  inline double Zq(double x) { return 1.0 / (1.0 + 0.2 * x) + 0.85; }
  inline double Mq(double x) { return 0.3 / (1.0 + 0.1 * x) + 0.05; }
  inline double dtZq(double x) { return -0.18 / (1.0 + 0.25 * x); }
  // meson-sector dressings (quark-meson flows, e.g. the hSigL sigma self-energy): the sigma/pion
  // wave-function (zSigL/zPiL) and the sigma/pion-quark Yukawa (hSigL/hPiL). Smooth positive stubs.
  inline double hSigL(double x) { return 0.9 + 0.3 * std::cos(0.45 * x); }
  inline double hPiL(double x) { return 0.8 + 0.35 * std::sin(0.4 * x); }
  inline double zSigL(double x) { return 1.0 / (1.0 + 0.15 * x * x) + 0.6; }
  inline double zPiL(double x) { return 1.0 / (1.0 + 0.2 * x * x) + 0.7; }
} // namespace shimdress

// Bundle of interpolators. The first seven are the ORDER the gauge-sector generated
// kernels expect (ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA); the quark-sector dressings
// (ZAqbq1, Zq, Mq, dtZq) are added for the Zq kernel.
struct DressingSet {
  DiFfRG::Fn ZA3{&shimdress::ZA3};
  DiFfRG::Fn ZAcbc{&shimdress::ZAcbc};
  DiFfRG::Fn ZA4{&shimdress::ZA4};
  DiFfRG::Fn dtZc{&shimdress::dtZc};
  DiFfRG::Fn Zc{&shimdress::Zc};
  DiFfRG::Fn dtZA{&shimdress::dtZA};
  DiFfRG::Fn ZA{&shimdress::ZA};
  DiFfRG::Fn ZAqbq1{&shimdress::ZAqbq1};
  DiFfRG::Fn ZAqbq4{&shimdress::ZAqbq4};
  DiFfRG::Fn ZAqbq7{&shimdress::ZAqbq7};
  DiFfRG::Fn Zq{&shimdress::Zq};
  DiFfRG::Fn Mq{&shimdress::Mq};
  DiFfRG::Fn dtZq{&shimdress::dtZq};
  DiFfRG::Fn hSigL{&shimdress::hSigL};
  DiFfRG::Fn hPiL{&shimdress::hPiL};
  DiFfRG::Fn zSigL{&shimdress::zSigL};
  DiFfRG::Fn zPiL{&shimdress::zPiL};
};
