#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/dense/dtensor.hpp"

namespace numtracer_kernels
{
  template<typename REG>
  class Sigl_dirac_dense
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double renv[3] = {};
      renv[0] = cos1 * l1;
      renv[1] = sqrt(1. - powr<2>(cos1)) * l1;
      renv[2] = p;
      const double _tr0_0 = numtracer::dense::contract_all(numtracer::dense::identity<100, 101>(), numtracer::dense::gamma_axis<1, 104, 101>(), numtracer::dense::gamma5<103, 104>(), numtracer::dense::gamma_axis<0, 102, 103>(), numtracer::dense::gamma5<100, 102>(), numtracer::dense::vector<0, 0, 3>(renv), numtracer::dense::vector<1, 2, 1>(renv)).scalar_value().re;
      const auto _interp1 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RF(powr<2>(k), powr<2>(p));
      const auto _den1 = powr<-1>(_interp2 + powr<2>(p));
      const auto _den2 = powr<-1>(_interp1 + powr<2>(l1));
      return -0.25 * _den1 * _den2 * _tr0_0;
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
using numtracer_kernels::Sigl_dirac_dense;