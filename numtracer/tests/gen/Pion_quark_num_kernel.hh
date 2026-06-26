#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "Pion_quark_num_kernels.hh"

namespace numtracer_kernels
{
  template<typename REG>
  class Pion_quark_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[numtracer_kernels::pion_quark_num::nenv];
      numtracer_kernels::pion_quark_num::fill(fenv, l1, cos1, p);
      const auto _interp1 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp3 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _den1 = powr<-2>(_interp1 + powr<2>(l1));
      const auto _den2 = powr<-1>(_interp2 + powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      return 0.6666666666666666 * _den1 * _den2 * _interp3 * numtracer_kernels::pion_quark_num::tr0(fenv);
    }

    static inline auto constant(const double& p, const double& k)
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
using numtracer_kernels::Pion_quark_num_kernel;