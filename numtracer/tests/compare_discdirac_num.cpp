// Codegen gate for a DISCONNECTED Dirac x Lorentz diagram: validate the NUMERIC backend kernel
// generated for a hand-built net whose diagrams each factorize into a closed Dirac trace, a constant
// colour trace, and an index-disjoint pure-Lorentz scalar (see tests/gen/gen_discdirac_numeric.wls)
// against the DENSE oracle, over random kinematics. This is the regression guard for the
// MakeNTKernel::disconnectmix case (the hSigL sigma-quark vertex, QCD.nb): the numeric backend now
// FACTORS the disconnected Lorentz scalar as its own trace and multiplies it in at assembly
// (coeff * tr_dirac * tr_colour * tr_lorentz), matching dense. The kernels carry no dressings, so REG
// is an unused stub.
#include "Discdirac_num_dense_kernel.hh"
#include "Discdirac_num_kernel.hh"

#include <cmath>
#include <cstdio>
#include <random>

struct Reg {}; // no regulator/dressing is referenced by these kernels

int main()
{
  using Dense = numtracer_kernels::Discdirac_num_dense_kernel<Reg>;
  using Num = numtracer_kernels::Discdirac_num_kernel<Reg>;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> Ul(0.05, 3.0), Uc(-0.999, 0.999);
  // This kernel is a SUM of products of independent traces (d1 + d2 + d3); the summands can nearly
  // cancel, so a pure pointwise-relative metric blows up at zeros of the sum even when numeric and
  // dense agree to machine precision. Gate on the max ABSOLUTE error relative to the kernel's overall
  // scale (maxAbs) — the standard, cancellation-robust check.
  double maxAbsErr = 0, maxAbs = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng), cos1 = Uc(rng), cos2 = Uc(rng), p = Ul(rng);
    const double d = Dense::kernel(l1, cos1, cos2, p);
    const double n = Num::kernel(l1, cos1, cos2, p);
    maxAbsErr = std::max(maxAbsErr, std::fabs(n - d));
    maxAbs = std::max(maxAbs, std::fabs(d));
  }
  const double rel = maxAbsErr / (1e-300 + maxAbs);
  std::printf("disconnected Dirac x Lorentz codegen numeric vs dense:  maxAbsErr=%.3e  maxAbs=%.3e  scaled=%.3e\n",
              maxAbsErr, maxAbs, rel);
  // maxAbs guards against a degenerate all-zero kernel masking a real mismatch.
  const bool ok = rel < 1e-12 && maxAbs > 1e-6;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
