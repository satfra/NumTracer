// NumTracer ZA-bench — standalone shim so the DiFfRG-targeted ZA kernels (both the
// FormTracer oracle and the NumTracer-generated kernel) compile and run with NO DiFfRG /
// Kokkos / CUDA dependency. It provides every name the generated kernels reference:
//   * the CUDA decorators (__host__/__device__) and Kokkos macro, neutralised for a CPU build,
//   * the math names (powr, pow, sqrt, fma) in DiFfRG::compute,
//   * the interpolator / coordinate / memory type stubs,
//   * a fixed analytic Regulator and the FULL set of 15 QCD dressing functions, used
//     IDENTICALLY by both evaluation paths.
//
// The specific regulator / dressing formulas are arbitrary but fixed: the benchmark checks
// that the two PATHS (FormTracer vs NumTracer) agree as algebraic functions, so any consistent
// definitions on both sides suffice. They reuse the formulas of numtracer/tests/refshim/shim.hpp
// so the two harnesses are comparable.
#pragma once

// The QCD kernels carry the GPU decorator "static __host__ __device__ inline"; neutralise the
// CUDA attributes for a host compile. (Kokkos kernels use KOKKOS_FORCEINLINE_FUNCTION instead.)
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#ifndef KOKKOS_FORCEINLINE_FUNCTION
#define KOKKOS_FORCEINLINE_FUNCTION inline
#endif

#include <cmath>
#include <complex>

namespace DiFfRG
{
  // Complex type the generated kernels resolve unqualified (emitted without `using std::complex;`
  // so a CUDA build can substitute cuda::std::complex). Here it is plain std::complex.
  template <class T> using complex = std::complex<T>;

  namespace compute
  {
    // compile-time integer power, matching DiFfRG::compute::powr.
    template <int n, typename NumberType>
      requires requires(NumberType x) { x *x; NumberType(1.) / x; }
    constexpr inline NumberType powr(const NumberType x)
    {
      if constexpr (n == 0) return NumberType(1.);
      else if constexpr (n < 0) return NumberType(1.) / powr<-n, NumberType>(x);
      else if constexpr (n == 1) return x;
      else if constexpr (n % 2 == 0) return powr<n / 2>(x) * powr<n / 2>(x);
      else return powr<n / 2>(x) * powr<n / 2>(x) * x;
    }
    using std::fma;
    using std::pow;
    using std::sqrt;
    inline std::complex<double> fma(const std::complex<double> &a, const std::complex<double> &b,
                                    const std::complex<double> &c)
    {
      return a * b + c;
    }
  } // namespace compute

  // type stubs the generated kernel signatures mention
  struct GPU_memory {};
  template <class T> struct LogarithmicCoordinates1D {};

  // a 1-D "interpolator" wrapping a function pointer (double -> double).
  struct Fn {
    double (*f)(double) = nullptr;
    double operator()(double x) const { return f(x); }
  };
  template <class A, class B, class C> using SplineInterpolator1D = Fn;

  // fixed analytic regulator (same on both paths; formulas from refshim/shim.hpp).
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
// The first 11 reuse refshim/shim.hpp verbatim; the full QCD ZA adds ZAqbq4/7 and ZAAqbq1/2.
namespace shimdress
{
  inline double ZA(double x) { return 1.0 / (1.0 + 0.5 * x * x) + 0.7; }
  inline double Zc(double x) { return 1.0 / (1.0 + 0.3 * x) + 0.9; }
  inline double ZAcbc(double x) { return 1.2 + 0.4 * std::sin(0.5 * x); }
  inline double ZA3(double x) { return 1.1 + 0.2 * std::cos(0.3 * x); }
  inline double ZA4(double x) { return 1.05 + 0.1 * std::sin(0.7 * x); }
  inline double dtZA(double x) { return -0.2 / (1.0 + x * x); }
  inline double dtZc(double x) { return -0.15 / (1.0 + 0.2 * x); }
  inline double ZAqbq1(double x) { return 1.15 + 0.25 * std::cos(0.4 * x); }
  inline double Zq(double x) { return 1.0 / (1.0 + 0.2 * x) + 0.85; }
  inline double Mq(double x) { return 0.3 / (1.0 + 0.1 * x) + 0.05; }
  inline double dtZq(double x) { return -0.18 / (1.0 + 0.25 * x); }
  // extra quark-gluon / two-gluon-two-quark vertex dressings of the full AAqbq ZA flow
  inline double ZAqbq4(double x) { return 1.08 + 0.18 * std::sin(0.35 * x); }
  inline double ZAqbq7(double x) { return 0.95 + 0.22 * std::cos(0.45 * x); }
  inline double ZAAqbq1(double x) { return 1.12 + 0.15 * std::cos(0.5 * x); }
  inline double ZAAqbq2(double x) { return 0.9 + 0.2 * std::sin(0.6 * x); }
} // namespace shimdress

// 15-dressing bundle in the ORDER the full-QCD ZA kernels expect:
//   ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, ZAAqbq1, ZAAqbq2, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq
struct DressingSetZA {
  DiFfRG::Fn ZA3{&shimdress::ZA3};
  DiFfRG::Fn ZAcbc{&shimdress::ZAcbc};
  DiFfRG::Fn ZA4{&shimdress::ZA4};
  DiFfRG::Fn ZAqbq1{&shimdress::ZAqbq1};
  DiFfRG::Fn ZAqbq4{&shimdress::ZAqbq4};
  DiFfRG::Fn ZAqbq7{&shimdress::ZAqbq7};
  DiFfRG::Fn ZAAqbq1{&shimdress::ZAAqbq1};
  DiFfRG::Fn ZAAqbq2{&shimdress::ZAAqbq2};
  DiFfRG::Fn dtZc{&shimdress::dtZc};
  DiFfRG::Fn Zc{&shimdress::Zc};
  DiFfRG::Fn dtZA{&shimdress::dtZA};
  DiFfRG::Fn ZA{&shimdress::ZA};
  DiFfRG::Fn dtZq{&shimdress::dtZq};
  DiFfRG::Fn Zq{&shimdress::Zq};
  DiFfRG::Fn Mq{&shimdress::Mq};
};
