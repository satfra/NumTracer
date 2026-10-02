// step-09 — Every codegen variant computes the same number.
//
// `options.wls` emitted the SAME network three ways. Almost every MakeNTKernel option is meant to be
// semantics-preserving: it changes how the kernel is spelled, not what it returns. This program is
// the assertion of that claim.
//
// It also checks the one option that DOES change observable behaviour — "Constant", which populates
// the separate constant() entry point — and confirms it leaves kernel() alone.
#include "Opt_const_kernel.hh"
#include "Opt_default_kernel.hh"
#include "Opt_ns_kernel.hh"

#include <cmath>
#include <cstdio>

namespace {
int failures = 0;

// Any callable works: the kernel takes the dressing as `const auto&` and only ever calls it.
struct MyZ {
  double operator()(double x) const { return 1.0 + 0.5 * x * x; }
};

void same(const char *what, double got, double ref) {
  const bool ok = std::fabs(got - ref) <= 1e-13 * std::max(1.0, std::fabs(ref));
  if (!ok) ++failures;
  std::printf("  %-28s = %14.10g   %s\n", what, got, ok ? "== default" : "<-- DIFFERS");
}
} // namespace

int main() {
  const MyZ myZ;
  const double p = 1.7, l1 = 0.9, cos1 = 0.35;

  // Closed form: the P^T and P^L diagrams add to the plain metric case (q1.q1 = p^2), both carrying
  // myZ(l1); the third diagram is (q1.ql) * q1.q1.
  const double want = myZ(l1) * p * p + (p * l1 * cos1) * p * p;
  const double ref = Opt_default_kernel::kernel(l1, cos1, p, myZ);

  std::printf("one network, three emissions, at (p, l1, cos1) = (%g, %g, %g)\n", p, l1, cos1);
  std::printf("  %-28s = %14.10g   (closed form %.10g)\n", "default", ref, want);
  if (std::fabs(ref - want) > 1e-12 * std::max(1.0, std::fabs(want))) {
    ++failures;
    std::printf("  default kernel disagrees with the closed form!\n");
  }

  same("Constant -> myZ[p]", Opt_const_kernel::kernel(l1, cos1, p, myZ), ref);
  same("renamed namespaces", Opt_ns_kernel::kernel(l1, cos1, p, myZ), ref);

  // "Constant" is the one option above that changes observable behaviour, and it changes only
  // constant() — the loop-independent term the integrator flat-adds after quadrature.
  std::printf("the constant() entry point\n");
  std::printf("  %-28s = %14.10g   (expected 0)\n", "default constant(p)",
              Opt_default_kernel::constant(p, myZ));
  std::printf("  %-28s = %14.10g   (expected myZ(p) = %.10g)\n", "Constant -> myZ[p]",
              Opt_const_kernel::constant(p, myZ), myZ(p));
  if (std::fabs(Opt_default_kernel::constant(p, myZ)) > 1e-15) ++failures;
  if (std::fabs(Opt_const_kernel::constant(p, myZ) - myZ(p)) > 1e-13) ++failures;

  std::printf(failures == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return failures == 0 ? 0 : 1;
}
