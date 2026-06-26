#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/dense/dtensor.hpp"

namespace numtracer_kernels
{
  template<typename REG>
  class Proj_num_dense_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double renv[7] = {};
      renv[0] = cos1 * l1;
      renv[1] = sqrt(1. - powr<2>(cos1)) * cos2 * l1;
      renv[2] = sqrt((-1. + powr<2>(cos1)) * (-1. + powr<2>(cos2))) * l1;
      renv[3] = p;
      renv[4] = -0.5 * p;
      renv[5] = 0.866025403784439 * p;
      renv[6] = 1.0 / (powr<2>(l1));
      const double _tr0_0 = numtracer::dense::contract_all(numtracer::dense::longitudinal_projector<0, 1, 0, 7, 6>(renv), numtracer::dense::vector<0, 3, 1>(renv), numtracer::dense::vector<1, 3, 1>(renv)).scalar_value().re;
      const double _tr1_0 = numtracer::dense::contract_all(numtracer::dense::transverse_projector<0, 1, 0, 7, 6>(renv), numtracer::dense::vector<0, 4, 3>(renv), numtracer::dense::vector<1, 4, 3>(renv)).scalar_value().re;  return _tr0_0 + _tr1_0;
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
using numtracer_kernels::Proj_num_dense_kernel;