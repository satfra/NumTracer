#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Cplxdiscf_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace numtracer_kernels
{
  class Cplxdiscf_num_kernel
  {
    public:
    #if NT_CPLXDISCF_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    static inline auto kernel(const float& l1, const float& cos1, const float& cos2, const float& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      float fenv[(numtracer_kernels::cplxdiscf_num::nenv) > 0 ? (numtracer_kernels::cplxdiscf_num::nenv) : 1];
      numtracer_kernels::cplxdiscf_num::fill(fenv, l1, cos1, cos2, p);
      return 0.f;
    }
    #elif NT_CPLXDISCF_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    static inline auto kernel(const float& l1, const float& cos1, const float& cos2, const float& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      float fenv[(numtracer_kernels::cplxdiscf_num::nenv) > 0 ? (numtracer_kernels::cplxdiscf_num::nenv) : 1];
      numtracer_kernels::cplxdiscf_num::fill(fenv, l1, cos1, cos2, p);
      const auto _interp1 = ntIm(numtracer_kernels::cplxdiscf_num::tr0(fenv));
      const auto _interp2 = ntIm(numtracer_kernels::cplxdiscf_num::tr1(fenv));
      const auto _interp3 = ntIm(numtracer_kernels::cplxdiscf_num::tr3(fenv));
      const auto _interp4 = ntIm(numtracer_kernels::cplxdiscf_num::tr4(fenv));
      const auto _interp5 = ntIm(numtracer_kernels::cplxdiscf_num::tr5(fenv));
      const auto _interp6 = ntIm(numtracer_kernels::cplxdiscf_num::tr6(fenv));
      const auto _interp7 = ntIm(numtracer_kernels::cplxdiscf_num::tr7(fenv));
      const auto _interp8 = ntRe(numtracer_kernels::cplxdiscf_num::tr2(fenv));
      const auto _interp9 = ntRe(numtracer_kernels::cplxdiscf_num::tr1(fenv));
      const auto _interp10 = ntRe(numtracer_kernels::cplxdiscf_num::tr3(fenv));
      const auto _interp11 = ntRe(numtracer_kernels::cplxdiscf_num::tr4(fenv));
      const auto _interp12 = ntIm(numtracer_kernels::cplxdiscf_num::tr2(fenv));
      const auto _interp13 = ntRe(numtracer_kernels::cplxdiscf_num::tr5(fenv));
      const auto _interp14 = ntRe(numtracer_kernels::cplxdiscf_num::tr6(fenv));
      const auto _interp15 = ntRe(numtracer_kernels::cplxdiscf_num::tr7(fenv));
      return fma(-1.f, _interp1, fma(-1.f, _interp12 * _interp13 * _interp14 * _interp15, fma(-1.f, _interp10 * _interp11 * _interp2, fma(_interp2, _interp3 * _interp4, fma(_interp12, _interp15 * _interp5 * _interp6, fma(_interp12, _interp14 * _interp5 * _interp7, fma(_interp12, _interp13 * _interp6 * _interp7, fma(-1.f, _interp14 * _interp15 * _interp5 * _interp8, fma(-1.f, _interp13 * _interp15 * _interp6 * _interp8, fma(-1.f, _interp13 * _interp14 * _interp7 * _interp8, fma(_interp5, _interp6 * _interp7 * _interp8, fma(-1.f, _interp11 * _interp3 * _interp9, fma(-1.f, _interp10 * _interp4 * _interp9, 0.f)))))))))))));
    }
    #else                              // the imaginary part survives: genuinely complex
    static inline auto kernel(const float& l1, const float& cos1, const float& cos2, const float& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      float fenv[(numtracer_kernels::cplxdiscf_num::nenv) > 0 ? (numtracer_kernels::cplxdiscf_num::nenv) : 1];
      numtracer_kernels::cplxdiscf_num::fill(fenv, l1, cos1, cos2, p);
      return complex<float>(0.f,1.f) * fma(numtracer_kernels::cplxdiscf_num::tr1(fenv), numtracer_kernels::cplxdiscf_num::tr3(fenv) * numtracer_kernels::cplxdiscf_num::tr4(fenv), fma(numtracer_kernels::cplxdiscf_num::tr2(fenv), numtracer_kernels::cplxdiscf_num::tr5(fenv) * numtracer_kernels::cplxdiscf_num::tr6(fenv) * numtracer_kernels::cplxdiscf_num::tr7(fenv), numtracer_kernels::cplxdiscf_num::tr0(fenv)));
    }
    #endif

    static inline auto constant(const float& p)
    {
      return 0.f;
    }
    private:
    static inline float ntRe(float x) { return x; }
    static inline float ntIm(float) { return 0.f; }
    template <class T> static inline auto ntRe(const T &z) -> decltype(z.real()) { return z.real(); }
    template <class T> static inline auto ntRe(const T &z) -> decltype(real(z)) requires (!requires { z.real(); }) { return real(z); }
    template <class T> static inline auto ntIm(const T &z) -> decltype(z.imag()) { return z.imag(); }
    template <class T> static inline auto ntIm(const T &z) -> decltype(imag(z)) requires (!requires { z.imag(); }) { return imag(z); }
  };
}
using numtracer_kernels::Cplxdiscf_num_kernel;