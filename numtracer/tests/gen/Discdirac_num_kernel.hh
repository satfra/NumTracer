#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "Discdirac_num_kernels.hh"

namespace numtracer_kernels
{
  template<typename REG>
  class Discdirac_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[numtracer_kernels::discdirac_num::nenv];
      numtracer_kernels::discdirac_num::fill(fenv, l1, cos1, cos2, p);  return fma(numtracer_kernels::discdirac_num::tr0(fenv), numtracer_kernels::discdirac_num::tr3(fenv) * numtracer_kernels::discdirac_num::tr4(fenv), fma(numtracer_kernels::discdirac_num::tr1(fenv), numtracer_kernels::discdirac_num::tr5(fenv) * numtracer_kernels::discdirac_num::tr6(fenv), fma(numtracer_kernels::discdirac_num::tr2(fenv), numtracer_kernels::discdirac_num::tr7(fenv) * numtracer_kernels::discdirac_num::tr8(fenv), 0.)));
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }private: static inline auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
    static inline auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
    static inline auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
    static inline auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
    static inline auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
    static inline auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }
  };
}
using numtracer_kernels::Discdirac_num_kernel;