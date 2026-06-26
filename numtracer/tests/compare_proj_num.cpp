// Codegen gate for the longitudinal projector: validate the NUMERIC backend kernel generated for a
// hand-built net mixing ntLongProj + ntTransProj (see tests/gen/gen_proj_numeric.wls) against the
// DENSE oracle (which builds longitudinal_projector / transverse_projector directly), over random
// kinematics. The kernels carry no dressings, so REG is an unused stub.
#include "Proj_num_dense_kernel.hh"
#include "Proj_num_kernel.hh"

#include <cmath>
#include <cstdio>
#include <random>

struct Reg {}; // no regulator/dressing is referenced by these kernels

int main()
{
  using Dense = numtracer_kernels::Proj_num_dense_kernel<Reg>;
  using Num = numtracer_kernels::Proj_num_kernel<Reg>;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> Ul(0.05, 3.0), Uc(-0.999, 0.999);
  double pw = 0, maxAbs = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng), cos1 = Uc(rng), cos2 = Uc(rng), p = Ul(rng);
    const double d = Dense::kernel(l1, cos1, cos2, p);
    const double n = Num::kernel(l1, cos1, cos2, p);
    pw = std::max(pw, std::fabs(n - d) / (1e-300 + std::fabs(d)));
    maxAbs = std::max(maxAbs, std::fabs(d));
  }
  std::printf("projector codegen numeric vs dense:  pointwise=%.3e  maxAbs=%.3e\n", pw, maxAbs);
  // maxAbs guards against a degenerate all-zero kernel masking a real mismatch.
  const bool ok = pw < 1e-10 && maxAbs > 1e-6;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
