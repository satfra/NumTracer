// NumTracer — on-the-fly color-contraction helper (the Wolfram<->C++ seam).
//
// In the full pipeline, Wolfram emits a tiny program like this for each distinct
// color network, compiles & runs it, and reads the printed number back in
// (Wolfram's symbolic color algebra is too slow). Here it computes the ghost-loop
// color factor f^{a c d} f^{b c d} for SU(N) and prints it.
//
// N is a runtime argument: the SU(N) data comes from the same runtime source the engine uses
// (sun_data_for: the typed tables for N = 2, 3, the generalized-Gell-Mann construction otherwise).
//
// Usage:  sun_contract <N>           # prints the diagonal value (= C_A = N)
//         sun_contract <N> --full    # prints the full a,b matrix (diagnostic)
#include "numtracer/network/sun_net.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int argc, char **argv) {
  const int N = (argc > 1) ? std::atoi(argv[1]) : 3;
  const bool full = (argc > 2) && std::strcmp(argv[2], "--full") == 0;
  if (N < 2) {
    std::fprintf(stderr, "N must be >= 2\n");
    return 1;
  }
  const auto &f = numtracer::network::sun_net_detail::sun_data_for(N).f_nz;
  const int A = N * N - 1;
  // f^{acd} f^{bcd} contracted over (c,d) — an inline numeric contraction (the library folds
  // colour numerically via sun_value_cx).
  auto ff = [&](int a, int b) {
    double s = 0;
    for (const auto &e1 : f)
      if (e1.a == a)
        for (const auto &e2 : f)
          if (e2.a == b && e1.b == e2.b && e1.c == e2.c) s += e1.v * e2.v;
    return s;
  };
  if (full) {
    std::printf("# f^{acd} f^{bcd} for SU(%d), adjoint dim %d\n", N, A);
    for (int a = 0; a < A; ++a) {
      for (int b = 0; b < A; ++b)
        std::printf("%8.4f ", ff(a, b));
      std::printf("\n");
    }
    return 0;
  }
  double maxoff = 0, mindiag = 1e300, maxdiag = -1e300;
  for (int a = 0; a < A; ++a)
    for (int b = 0; b < A; ++b) {
      const double v = ff(a, b);
      if (a == b) {
        mindiag = std::min(mindiag, v);
        maxdiag = std::max(maxdiag, v);
      } else {
        maxoff = std::max(maxoff, std::fabs(v));
      }
    }
  if (maxoff > 1e-9 || (maxdiag - mindiag) > 1e-9) {
    std::fprintf(stderr, "non-scalar result: off=%.3e diag=[%.6f,%.6f]\n", maxoff, mindiag, maxdiag);
    return 2;
  }
  std::printf("%.15g\n", maxdiag); // the single exact number Wolfram reads back
  return 0;
}
