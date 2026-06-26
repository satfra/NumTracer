// CPU same-quadrature reference for the ZA4 GPU integration test (see ref.hpp / ref_za3.cpp).
#include "ZA4_num_kernel.hh"
#include "ref.hpp"

#include <complex>

std::vector<double> cpu_reference_za4(const Axis ax[4], const std::vector<double> &p_grid, double k,
                                      double prefactor) {
  using K = DiFfRG::ZA4_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  std::vector<double> out(p_grid.size());
  for (size_t ip = 0; ip < p_grid.size(); ++ip) {
    const double p = p_grid[ip];
    long double acc = 0.0L;
    for (size_t i = 0; i < ax[0].nodes.size(); ++i) {
      const double l1 = ax[0].nodes[i];
      for (size_t j = 0; j < ax[1].nodes.size(); ++j)
        for (size_t l = 0; l < ax[2].nodes.size(); ++l)
          for (size_t m = 0; m < ax[3].nodes.size(); ++m) {
            const double w = prefactor * ax[0].weights[i] * ax[1].weights[j] * ax[2].weights[l] *
                             ax[3].weights[m];
            const double v = std::real(K::kernel(l1, ax[1].nodes[j], ax[2].nodes[l], ax[3].nodes[m], p, k,
                                                 d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA,
                                                 d.ZA, d.dtZq, d.Zq, d.Mq));
            acc += (long double)(w * l1 * l1 * l1 * v);
          }
    }
    acc += std::real(K::constant(p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA,
                                 d.dtZq, d.Zq, d.Mq));
    out[ip] = (double)acc;
  }
  return out;
}
