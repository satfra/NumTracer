// Regression gate for the PACKED renv layout of the dense momentum builders.
//
// The dense builders store only a momentum's structurally-nonzero components, packed: the k-th
// present component (in ascending Lorentz order) lives at renv[Base+k]. The committed kernels and
// tutorials only exercise contiguous-from-bit-0 masks (1, 3, 7), for which packed indexing equals
// the old index-by-Lorentz-component scheme. This test exercises the case they don't: a GAPPED
// mask (0b1001 — components 0 and 3 present, 1 and 2 absent), where Base+k != Base+component and
// the packing actually matters.
#include "numtracer/dense/dtensor.hpp"

#include <cmath>
#include <cstdio>

using namespace numtracer;

enum { mu, nu };

int main() {
  int fails = 0;

  // q = (2.5, 0, 0, -1.5): mask 0b1001, the two present components packed into renv[0], renv[1].
  const double q0 = 2.5, q3 = -1.5;
  double renv[4] = {};
  renv[0] = q0; // -> Lorentz component 0
  renv[1] = q3; // -> Lorentz component 3 (packed: 2nd present component)
  renv[2] = 1.0 / (q0 * q0 + q3 * q3); // 1/q^2 slot

  // (1) vector: packed slots must land on Lorentz components 0 and 3, with 1 and 2 zero.
  {
    auto v = dense::vector<mu, 0, 0b1001>(renv);
    const double c[4] = {v.data[0].re, v.data[1].re, v.data[2].re, v.data[3].re};
    const bool ok = c[0] == q0 && c[1] == 0.0 && c[2] == 0.0 && c[3] == q3;
    std::printf("  vector<0b1001>  (%.2f, %.2f, %.2f, %.2f)  %s\n", c[0], c[1], c[2], c[3], ok ? "ok" : "FAIL");
    if (!ok) ++fails;
  }

  // (2) transverse projector P_{ij}(q) = delta_ij - q_i q_j / q^2 with the same gapped mask.
  {
    auto P = dense::transverse_projector<mu, nu, 0, 0b1001, 2>(renv);
    const double iq2 = renv[2];
    double q[4] = {q0, 0.0, 0.0, q3};
    double err = 0.0;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const double want = (i == j ? 1.0 : 0.0) - q[i] * q[j] * iq2;
        err = std::max(err, std::fabs(P.data[static_cast<std::size_t>(i) * 4 + j].re - want));
      }
    const bool ok = err < 1e-14;
    std::printf("  transverse_projector<0b1001>  max|err|=%.2e  %s\n", err, ok ? "ok" : "FAIL");
    if (!ok) ++fails;
  }

  // (3) longitudinal projector P^L_{ij}(q) = q_i q_j / q^2, gapped mask.
  {
    auto P = dense::longitudinal_projector<mu, nu, 0, 0b1001, 2>(renv);
    const double iq2 = renv[2];
    double q[4] = {q0, 0.0, 0.0, q3};
    double err = 0.0;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        err = std::max(err, std::fabs(P.data[static_cast<std::size_t>(i) * 4 + j].re - q[i] * q[j] * iq2));
    const bool ok = err < 1e-14;
    std::printf("  longitudinal_projector<0b1001>  max|err|=%.2e  %s\n", err, ok ? "ok" : "FAIL");
    if (!ok) ++fails;
  }

  std::printf(fails == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return fails ? 1 : 0;
}
