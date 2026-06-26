// CPU same-quadrature reference for the ZA3 GPU integration test (see ref.hpp). Uses the
// committed HOST kernel + the shared tests/gpu shim (this target's include path puts
// tests/gpu first, so the kernel's `#include "shim.hpp"` resolves to the host/device shim —
// identical formulas to tests/refshim).
#include "ZA3_num_kernel.hh"
#include "ref.hpp"

#include <complex>

std::vector<double> cpu_reference_za3(const Axis ax[3], const std::vector<double> &p_grid, double k,
                                      double prefactor) {
  using K = DiFfRG::ZA3_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  std::vector<double> out(p_grid.size());
  for (size_t ip = 0; ip < p_grid.size(); ++ip) {
    const double p = p_grid[ip];
    long double acc = 0.0L;
    for (size_t i = 0; i < ax[0].nodes.size(); ++i) {
      const double l1 = ax[0].nodes[i];
      for (size_t j = 0; j < ax[1].nodes.size(); ++j)
        for (size_t l = 0; l < ax[2].nodes.size(); ++l) {
          const double w = prefactor * ax[0].weights[i] * ax[1].weights[j] * ax[2].weights[l];
          const double v = std::real(K::kernel(l1, ax[1].nodes[j], ax[2].nodes[l], p, k, d.ZA3, d.ZAcbc,
                                               d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq,
                                               d.Mq));
          acc += (long double)(w * l1 * l1 * l1 * v);
        }
    }
    acc += std::real(K::constant(p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA,
                                 d.dtZq, d.Zq, d.Mq));
    out[ip] = (double)acc;
  }
  return out;
}
