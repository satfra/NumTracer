// NumTracer — ZAqbq1 (quark-gluon vertex) validation + timing.
//
// A 2-loop-angle Dirac VERTEX flow whose projector carries an `i`, so the physical flow is the
// real part of a complex integrand. TWO imaginary units combine to make it real: the projector's
// `i` and the non-abelian colour factor f^{abc}T^bT^c = (iN/2)T^a. The NUMERIC kernel keeps the
// heavy momentum trace real and folds the per-diagram colour×coefficient complex, returns
// std::complex<double>, and the caller takes std::real. (The earlier dense baseline overflowed the
// stack at runtime and the symbolic et kernel OOMs at compile time — the numeric matrix-product
// backend is the single engine path here.) Validated against the copied FormTracer oracle.
// Machine-readable RESULT lines feed the benchmark report.
#include "ZAqbq1_num_kernel.hh" // generated numeric kernel  (gen/)
#include "ZAqbq1_kernel.hh"     // copied FormTracer oracle  (refshim/)
#include "shim.hpp"

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <vector>

int main() {
  using Form = DiFfRG::ZAqbq1_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZAqbq1_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;
  const double Nf = 2.0;

  auto callForm = [&](double l1, double c1, double c2, double p, double k) {
    return Form::kernel(l1, c1, c2, p, k, Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq);
  };
  // the numeric kernel returns std::complex<double> (the vertex projector's i + imaginary
  // non-abelian colour); the physical flow is the real part. The heavy trace polynomials stay
  // real — only the per-diagram colour×coefficient is complex.
  auto callNum = [&](double l1, double c1, double c2, double p, double k) {
    return std::real(
        Num::kernel(l1, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);

  // ZAqbq1 matches the oracle POINTWISE (no odd-cos loop-routing term), like the Zq 2-point.
  const int Nt = 200000;
  std::vector<double> TL(Nt), TC1(Nt), TC2(Nt), TP(Nt), TK(Nt);
  for (int i = 0; i < Nt; ++i) { TL[i] = U(rng); TC1[i] = Uc(rng); TC2[i] = Uc(rng); TP[i] = U(rng); TK[i] = U(rng); }

  double e_num = 0;
  for (int i = 0; i < Nt; ++i) {
    const double f = callForm(TL[i], TC1[i], TC2[i], TP[i], TK[i]);
    e_num = std::max(e_num, std::fabs(callNum(TL[i], TC1[i], TC2[i], TP[i], TK[i]) - f) / (1e-300 + std::fabs(f)));
  }
  std::printf("ZAqbq1 (quark-gluon vertex) vs Form, pointwise over %d points:\n", Nt);
  std::printf("  numeric max rel err = %.3e\n", e_num);

  auto bench = [&](auto fn, int n) {
    volatile double sink = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; ++i) sink += fn(i);
    auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / n;
  };
  const double t_form = bench([&](int i) { return callForm(TL[i], TC1[i], TC2[i], TP[i], TK[i]); }, Nt);
  const double t_num = bench([&](int i) { return callNum(TL[i], TC1[i], TC2[i], TP[i], TK[i]); }, Nt);
  std::printf("\ntiming:\n  Form    %12.2f ns/eval\n  numeric %12.2f ns/eval\n", t_form, t_num);
  std::printf("RESULT kernel=ZAqbq1 path=Form    ns=%.2f relerr=0\n", t_form);
  std::printf("RESULT kernel=ZAqbq1 path=numeric ns=%.2f relerr=%.3e\n", t_num, e_num);

  // 1e-8 is the numeric-frame tolerance used across the numeric kernels (compare_za4_num,
  // compare_za3_147_num): the fixed kinematic frame + the to_genprog noise-prune leave a ~1e-9
  // round-off floor, amplified here by the real-part cancellation of a complex integrand.
  const bool pass = e_num < 1e-8;
  std::printf("\nZAqbq1 (numeric vs Form): %s\n", pass ? "MATCH (< 1e-8)" : "MISMATCH");
  return pass ? 0 : 1;
}
