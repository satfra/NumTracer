// The typed-out SU(2)/SU(3) tables (sun/sun_data.hpp) against the generalized-Gell-Mann construction
// they were generated from (sun_net_detail::build_oracle, which also serves every other N). Production
// colour factors for N = 2, 3 come from the tables, so a wrong table entry would silently change
// every colour value; this test is what would notice.
//
// Bit-exact, with one allowance: the construction writes the -i/2 generator entries as -I*0.5, whose
// real part is -0.0, while the tables spell it +0.0. The two compare equal, and the colour fold
// treats them identically.
#include "numtracer/network/sun_net.hpp"

#include <cstdio>
#include <cstring>

namespace {
int fails = 0;

bool bitEqual(double x, double y) { return std::memcmp(&x, &y, sizeof(double)) == 0; }
// bit-equal, or both zero (of either sign)
bool sameValue(double x, double y) { return bitEqual(x, y) || (x == 0.0 && y == 0.0); }

template <int N> void check() {
  namespace snd = numtracer::network::sun_net_detail;
  const snd::SUNDyn tab = snd::seed_from_table<N>();
  const snd::SUNDyn ora = snd::build_oracle(N);
  const int adj = N * N - 1;
  int genDiff = 0, signedZeros = 0;
  for (int a = 0; a < adj; ++a)
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) {
        const auto t = tab.gens[a](i, j), o = ora.gens[a](i, j);
        if (!sameValue(t.real(), o.real()) || !sameValue(t.imag(), o.imag())) ++genDiff;
        else if (!bitEqual(t.real(), o.real()) || !bitEqual(t.imag(), o.imag())) ++signedZeros;
      }
  bool fOk = tab.f_nz.size() == ora.f_nz.size();
  for (std::size_t k = 0; fOk && k < tab.f_nz.size(); ++k) {
    const auto &x = tab.f_nz[k], &y = ora.f_nz[k];
    fOk = x.a == y.a && x.b == y.b && x.c == y.c && bitEqual(x.v, y.v);
  }
  std::printf("  SU(%d): generators %s (%d entries differ only in the sign of zero), f^{abc} %zu entries %s\n", N,
              genDiff == 0 ? "ok" : "FAIL", signedZeros, tab.f_nz.size(), fOk ? "bit-identical" : "FAIL");
  if (genDiff != 0 || !fOk) ++fails;
}
} // namespace

int main() {
  std::printf("typed SU(N) tables vs the generalized-Gell-Mann construction\n");
  check<2>();
  check<3>();
  std::printf(fails == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return fails == 0 ? 0 : 1;
}
