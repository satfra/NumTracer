#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "numtracer/expr/real_cse.hpp"
#include "ZA3_gpu_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZA3_gpu_kernel
  {
    public:
    // clang-format off
static __host__ __device__ inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
using namespace DiFfRG;
    // clang-format on

    using namespace DiFfRG::compute;
    using namespace numtracer;
    const double cosl1p1 = cos1;
    const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
    const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
    double fenv[DiFfRG::za3_gpu::nenv];
    const double dr_0 = Mq(l1);
    const double dr_1 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
    const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const double dr_4 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    DiFfRG::za3_gpu::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5);
    const auto _interp1 = ntRe(DiFfRG::za3_gpu::tr0(fenv));
    const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
    const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
    const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
    const auto _interp7 = ZA(l1);
    const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
    const auto _interp9 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp10 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
    const auto _interp11 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp12 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp13 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp14 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp15 = ntRe(DiFfRG::za3_gpu::tr1(fenv));
    const auto _interp16 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
    const auto _interp17 = ntRe(DiFfRG::za3_gpu::tr2(fenv));
    const auto _interp18 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp19 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp20 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp21 = Zc(k);
    const auto _interp22 = dtZc(k);
    const auto _interp23 = Zc(1.02 * k);
    const auto _interp24 = Zc(l1);
    const auto _interp25 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp26 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp27 = ntIm(DiFfRG::za3_gpu::tr3(fenv));
    const auto _interp28 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp29 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp30 = ZAqbq1(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp31 = RFdot(powr<2>(k), powr<2>(l1));
    const auto _interp32 = Zq(k);
    const auto _interp33 = RF(powr<2>(k), powr<2>(l1));
    const auto _interp34 = dtZq(k);
    const auto _interp35 = Zq(1.02 * k);
    const auto _interp36 = Mq(l1);
    const auto _interp37 = Zq(l1);
    const auto _interp38 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp39 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
    const auto _interp40 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _interp41 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _interp42 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
    const auto _interp43 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _den1 = powr<-1>(1. + powr<6>(k));
    const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp7 * powr<2>(l1));
    const auto _den3 = powr<-2>(_interp21 * _interp4 + _interp24 * powr<2>(l1));
    const auto _den4 = powr<-1>(_interp3 * _interp8 + _interp9 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den5 = -powr<-1>(_interp21 * _interp8 + _interp25 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den6 = powr<-2>(powr<2>(_interp36) + powr<2>(_interp32 * _interp33 + _interp37 * l1));
    const auto _den7 = powr<-1>(_interp10 * _interp3 + _interp11 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _den8 = -powr<-1>(_interp10 * _interp21 + _interp26 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _den9 = powr<-1>(powr<2>(_interp38) + powr<2>(_interp32 * _interp39 + _interp40 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
    const auto _den10 = powr<-1>(powr<2>(_interp41) + powr<2>(_interp43) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp32 * _interp42 * (_interp32 * _interp42 + 2. * _interp43 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
    const auto _cse1 = powr<-1>(l1);
    const auto _cse2 = powr<-2>(p);
    const auto _cse3 = -_interp3;
    const auto _cse4 = _cse3 + _interp6;
    const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
    const auto _cse6 = _cse5 + _interp5;
    const auto _cse7 = _cse6 * _interp4;
    const auto _cse8 = _interp2 * _interp3;
    const auto _cse9 = _cse7 + _cse8;
    return fma(-0.005050505050505051, _cse2 * _cse9 * _den2 * _den4 * _den7 * _interp1 * _interp12 * _interp13 * _interp14, fma(-0.005050505050505051, _cse2 * _cse9 * _den2 * _den7 * _interp13 * _interp15 * _interp16, fma(-0.0202020202020202, _cse2 * _den10 * _den6 * _den9 * _interp27 * _interp28 * _interp29 * _interp30 * (-_cse1 * _interp31 * _interp32 - _cse1 * _interp33 * (_interp34 + 50. * (-_interp32 + _interp35))), fma(0.0101010101010101, _cse2 * _den3 * _den5 * _den8 * _interp17 * _interp18 * _interp19 * _interp20 * (_interp2 * _interp21 + (_interp22 + 50. * (-_interp21 + _interp23)) * _interp4), 0.))));
  }

  // clang-format off
static __host__ __device__ inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
{
return 0.;
  // clang-format on

}private: static __host__ __device__ inline auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
static __host__ __device__ inline auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
static __host__ __device__ inline auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
static __host__ __device__ inline auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
static __host__ __device__ inline auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
static __host__ __device__ inline auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }

static __host__ __device__ inline double ntRe(double x) { return x; }
template <class T> static __host__ __device__ inline double ntRe(const T &z) { return z.real(); }
static __host__ __device__ inline double ntIm(double) { return 0.0; }
template <class T> static __host__ __device__ inline double ntIm(const T &z) { return z.imag(); }
};
}
using DiFfRG::ZA3_gpu_kernel;