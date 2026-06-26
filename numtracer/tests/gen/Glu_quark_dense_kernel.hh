#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/dense/dtensor.hpp"

namespace numtracer_kernels
{
  template<typename REG>
  class Glu_quark_dense
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double renv[4] = {};
      renv[0] = p;
      renv[1] = cos1 * l1;
      renv[2] = sqrt(1. - powr<2>(cos1)) * l1;
      renv[3] = 1.0 / (powr<2>(p));
      const double _tr0_0 = numtracer::dense::contract_all(numtracer::dense::gamma_axis<0, 100, 101>(), numtracer::dense::gamma_axis<3, 105, 100>(), numtracer::dense::gamma_axis<2, 104, 105>(), numtracer::dense::gamma_axis<4, 101, 102>(), numtracer::dense::gamma_axis<1, 102, 103>(), numtracer::dense::gamma_axis<5, 103, 104>(), numtracer::dense::transverse_projector<4, 5, 0, 1, 3>(renv), numtracer::dense::vector<0, 1, 3>(renv), numtracer::dense::vector<2, 1, 3>(renv), numtracer::dense::vector<3, 1, 3>(renv), numtracer::dense::add_all(numtracer::dense::vector<1, 1, 3>(renv), numtracer::dense::scale<numtracer::Cx{-1., 0.0}>(numtracer::dense::vector<1, 0, 1>(renv)))).scalar_value().re;
      constexpr double _tr0_1 = numtracer::dense::contract_all(numtracer::dense::delta_adj<3, 6, 7>(), numtracer::dense::T<3, 6, 9, 10>(), numtracer::dense::delta_fund<3, 8, 9>(), numtracer::dense::delta_fund<3, 10, 11>(), numtracer::dense::delta_fund<3, 13, 8>(), numtracer::dense::delta_fund<3, 12, 13>(), numtracer::dense::T<3, 7, 11, 12>()).scalar_value().re;
      const auto _interp1 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp3 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _den1 = powr<-2>(_interp1 + powr<2>(l1));
      const auto _den2 = powr<-1>(_interp2 + powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      return 0.1666666666666667 * _den1 * _den2 * _interp3 * _tr0_0 * _tr0_1;
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
using numtracer_kernels::Glu_quark_dense;