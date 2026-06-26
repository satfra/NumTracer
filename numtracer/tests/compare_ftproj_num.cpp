// Codegen gate for FINITE-TEMPERATURE support: validate the NUMERIC backend kernel generated for a
// finite-T net (see tests/gen/gen_ftproj_numeric.wls) — electric/magnetic/transverse projectors plus
// the spatial scalar product ntSPS and the temporal component ntVec[q,0], over the finite-T frame
// ftFrameQuark — against the DENSE oracle, over random kinematics. The temporal components p0 (e.g.
// the lowest fermionic Matsubara frequency πT) and l0 are independent runtime scalars here, so this
// exercises a genuinely broken-O(4) frame. The kernels carry no dressings, so REG is an unused stub.
#include "FTProj_num_dense_kernel.hh"
#include "FTProj_num_kernel.hh"

#include <cmath>
#include <cstdio>
#include <random>

struct Reg {}; // no regulator/dressing is referenced by these kernels

int main()
{
  using Dense = numtracer_kernels::FTProj_num_dense_kernel<Reg>;
  using Num = numtracer_kernels::FTProj_num_kernel<Reg>;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> Up(0.05, 3.0), Uc(-0.999, 0.999), Ut(-3.0, 3.0);
  double pw = 0, maxAbs = 0;
  for (int i = 0; i < 200000; ++i) {
    const double p0 = Ut(rng), p = Up(rng), l0 = Ut(rng), l1 = Up(rng), cos1 = Uc(rng);
    const double d = Dense::kernel(p0, p, l0, l1, cos1);
    const double n = Num::kernel(p0, p, l0, l1, cos1);
    pw = std::max(pw, std::fabs(n - d) / (1e-300 + std::fabs(d)));
    maxAbs = std::max(maxAbs, std::fabs(d));
  }
  std::printf("finite-T projector codegen numeric vs dense:  pointwise=%.3e  maxAbs=%.3e\n", pw, maxAbs);
  // maxAbs guards against a degenerate all-zero kernel masking a real mismatch.
  const bool ok = pw < 1e-10 && maxAbs > 1e-6;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
