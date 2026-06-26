// NumTracer GPU tests — quadrature rules from GSL (host-only; the CUDA integrator consumes
// the resulting node/weight arrays and never touches GSL itself).
//
// gsl_integration_fixed computes nodes/weights for ∫_a^b f(x) ω(x) dx ≈ Σ w_i f(x_i) with the
// rule's weight function ω baked into the w_i:
//   * legendre:   ω = 1            (radial momentum, cos2, phi axes)
//   * chebyshev2: ω = √((b-x)(x-a)) — on [-1,1] this is the √(1-cos²) angular jacobian of the
//                 4D measure, so the integrand must NOT multiply it (DiFfRG convention).
#pragma once

#include <gsl/gsl_integration.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

struct Axis {
  std::vector<double> nodes, weights;
};

enum class Quad { legendre, chebyshev2 };

inline Axis gsl_axis(Quad t, size_t order, double a, double b) {
  const gsl_integration_fixed_type *T =
      t == Quad::legendre ? gsl_integration_fixed_legendre : gsl_integration_fixed_chebyshev2;
  gsl_integration_fixed_workspace *w = gsl_integration_fixed_alloc(T, order, a, b, 0.0, 0.0);
  if (!w) {
    std::fprintf(stderr, "gsl_integration_fixed_alloc failed (order=%zu)\n", order);
    std::exit(2);
  }
  const double *n = gsl_integration_fixed_nodes(w);
  const double *wt = gsl_integration_fixed_weights(w);
  Axis ax{{n, n + order}, {wt, wt + order}};
  gsl_integration_fixed_free(w);
  return ax;
}
