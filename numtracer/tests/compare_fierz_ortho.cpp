// Projector-orthonormality gate for a FOUR-QUARK Fierz basis, through the full NumTracer pipeline.
// Kernel from tests/gen/gen_fierz_ortho_numeric.wls (basis FierzCompleteNf2Nc3NJL, n = 10).
//
// <P_i, V_j> = delta_ij over all 100 pairs, contracted in ONE net against a distinct-prime weight
// matrix w[i][j] = Prime[10(i-1)+j], so the kernel must equal
//     W = Sum_ij w_ij <P_i,V_j> = Sum_i w_ii = 2451.
// Distinct primes mean no plausible error cancels: one wrong pair shifts W by its own weight, and a
// uniform normalisation error shows up as an exact multiple (see below).
//
// This is a harder gate than the AqbqDirect ones: FOUR external legs and TWO independent spinor
// lines (so the multi-loop loop_sep path runs), a basis written in the finite-T 3+1 split (so every
// structure carries gamma^0, the fixed-Lorentz-component path), and structures 3,4,6,8,9,10 built
// from SU(N) Levi-Civita tensors in BOTH colour and flavour (the epsilon-pair contraction, including
// a pair straddling an eager Plus in the diquark vertex).
//
// Momentum-INDEPENDENT by construction: these are contact structures, so W must come out constant.
// The kernel is sampled over random momenta anyway — that constancy is itself part of the assertion.
//
// HISTORY, because the failure mode is instructive: this first read 9804 = 4 * 2451 exactly. The
// tracer was right; the GENERATOR SCRIPT was wrong. TBInfo[] reports this basis's inner product as
// "2 T1[1,2,3,4] T2[2,1,4,3] - 2 T1[1,2,3,4] T2[2,3,4,1]", but TBGetVertex already carries both leg
// orderings with opposite signs and an overall 1/2 — the vertex IS the antisymmetrised structure, so
// V_j(idxB) = -V_j(idxA) and forming that combination again gave 2F - 2(-F) = 4F. If this gate ever
// reads an exact small-integer multiple of 2451 again, suspect the net construction, not the engine.
#include "FierzOrtho_kernel.hh"

#include <cmath>
#include <cstdio>
#include <random>

int main()
{
  using K = numtracer_kernels::FierzOrtho_kernel;

  const double want = 2451.0;
  std::mt19937_64 rng(20260718);
  std::uniform_real_distribution<double> Ul(0.3, 2.5), Uc(-0.99, 0.99), Uphi(0.0, 6.283185307179586);

  double worst = 0.0, spread = 0.0, first = 0.0;
  const int N = 2000;
  for (int i = 0; i < N; ++i) {
    const double p = Ul(rng), l1 = Ul(rng), cos1 = Uc(rng), cos2 = Uc(rng), phi = Uphi(rng);
    const double v = K::kernel(p, l1, cos1, cos2, phi);
    if (i == 0) first = v;
    worst = std::max(worst, std::fabs(v - want));
    spread = std::max(spread, std::fabs(v - first));  // momentum independence, independent of `want`
  }

  std::printf("Fierz four-quark orthonormality (2451) vs analytic: max|W-2451|=%.3e  spread=%.3e\n",
              worst, spread);
  const bool ok = (worst < 1e-8) && (spread < 1e-8);
  if (!ok) {
    if (spread >= 1e-8)
      std::printf("TESTS FAILED (W depends on the momenta; these are CONTACT structures)\n");
    else
      std::printf("TESTS FAILED (W = %.10f, want 2451; an exact multiple points at the net "
                  "construction, not the engine)\n", first);
  } else {
    std::printf("ALL TESTS PASSED\n");
  }
  return ok ? 0 : 1;
}
