// step-06 — Consuming the kernel the front-end just generated.
//
// `first_kernel.wls` generated two headers into the build directory. This program is an ordinary
// consumer of them: it includes the kernel, calls it over a grid of kinematic points, and checks
// the answer against the closed form
//
//     p . P(l) . p / p^2  =  1 - cos^2(theta) .
//
// Note what this file does NOT do. It never mentions Frame, Poly, DiracChain or LorentzNet —
// none of the machinery of steps 01-05. A generated kernel is self-contained straight-line
// arithmetic that depends only on two small NumTracer headers (codegen/runtime.hpp and
// sun/sun_data.hpp); the engine ran at BUILD time and is not present at run time.
#include "first_kernel.hh" // the generated kernel class (pulls in first_kernels.hh + the runtime)

#include <cmath>
#include <cstdio>

int main() {
  double worst = 0.0;
  int bad = 0;

  // Sweep the two arguments the kernel actually depends on. `p` is swept too, to confirm it
  // cancels: the network carries p^2 from the two ntVec factors and 1/p^2 from the ntSP
  // coefficient, so the result must be p-independent.
  for (int ip = 1; ip <= 5; ++ip) {
    const double p = 0.25 * ip;
    for (int il = 1; il <= 5; ++il) {
      const double l1 = 0.4 * il;
      for (int ic = 0; ic <= 20; ++ic) {
        const double cos1 = -1.0 + 0.1 * ic;

        // The generated signature is exactly the "Args" list from the .wls, in order.
        const double got = first_kernel::kernel(l1, cos1, p);
        const double want = 1.0 - cos1 * cos1;

        const double err = std::fabs(got - want);
        worst = std::max(worst, err);
        if (err >= 1e-12) ++bad;
      }
    }
  }

  std::printf("generated kernel  p.P(l).p / p^2  vs  1 - cos^2(theta)\n");
  std::printf("  525 kinematic points, worst |error| = %.3e   (%d bad)\n", worst, bad);

  // The kernel has a second entry point, constant(): the loop-INDEPENDENT term that a flow adds
  // flat to the integral. This network has none, so the front-end emitted a literal 0.
  // step-13 is about the case where it is not zero.
  std::printf("  constant(p=1) = %g   (no loop-independent term in this network)\n",
              first_kernel::constant(1.0));

  const bool ok = (bad == 0);
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
