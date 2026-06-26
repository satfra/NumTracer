// PROTOTYPE (task #22): matrix-product contraction with multivariate-polynomial entries.
//
// The premise: γ/metric/projector are NUMERIC; only momenta (and dressings/denominators) are
// symbolic. So contract a diagram as 4×4 spinor matrix PRODUCTS whose entries are polynomials in the
// symbolic momentum data — bounded by the matrix structure, NOT the multilinear sp-monomial blowup
// that reduce_product suffers. This file validates the core machinery: a γ-trace built from
// polynomial-entry matrices, substituted to numbers, must equal a direct numeric 4×4 matrix trace.
//
// Build:  g++ -std=c++20 -O2 -I include tests/proto_matprod.cpp -o /tmp/proto && /tmp/proto
#include "numtracer/dirac/dirac_data.hpp" // kGamma (Euclidean Weyl, {g^mu,g^nu}=2 delta)
#include <array>
#include <cstdio>
#include <map>
#include <random>
#include <vector>

using numtracer::Cx;
using numtracer::dirac::kGamma;

// ---- multivariate polynomial over `nvars` symbolic variables ---------------------------------
// monomial = exponent vector (length nvars); poly = map<exponents, complex coeff>. Map-based: simple,
// auto-collects like terms (the whole point — collection happens continuously during the products).
struct MPoly {
  int nvars = 0;
  std::map<std::vector<int>, Cx> t; // exponents -> coeff
  MPoly() = default;
  explicit MPoly(int nv) : nvars(nv) {}
  static MPoly constant(int nv, Cx c) {
    MPoly p(nv);
    if (!(c.re == 0 && c.im == 0)) p.t[std::vector<int>(static_cast<std::size_t>(nv), 0)] = c;
    return p;
  }
  static MPoly var(int nv, int i) { // the i-th variable, coeff 1
    MPoly p(nv);
    std::vector<int> e(static_cast<std::size_t>(nv), 0);
    e[static_cast<std::size_t>(i)] = 1;
    p.t[e] = Cx{1, 0};
    return p;
  }
  std::size_t size() const { return t.size(); }
  void addTerm(const std::vector<int> &e, Cx c) {
    if (c.re == 0 && c.im == 0) return;
    auto it = t.find(e);
    if (it == t.end()) t[e] = c;
    else { it->second = it->second + c; if (it->second.re == 0 && it->second.im == 0) t.erase(it); }
  }
};
static MPoly operator+(const MPoly &a, const MPoly &b) {
  if (a.t.empty()) return b;
  if (b.t.empty()) return a;
  MPoly r = a;
  for (const auto &[e, c] : b.t) r.addTerm(e, c);
  return r;
}
static MPoly operator*(const MPoly &a, const MPoly &b) {
  if (a.t.empty() || b.t.empty()) return MPoly(a.nvars ? a.nvars : b.nvars);
  MPoly r(a.nvars);
  std::vector<int> e(static_cast<std::size_t>(a.nvars));
  for (const auto &[ea, ca] : a.t)
    for (const auto &[eb, cb] : b.t) {
      for (int k = 0; k < a.nvars; ++k) e[static_cast<std::size_t>(k)] = ea[static_cast<std::size_t>(k)] + eb[static_cast<std::size_t>(k)];
      r.addTerm(e, ca * cb);
    }
  return r;
}
static Cx eval(const MPoly &p, const std::vector<double> &x) {
  Cx s{0, 0};
  for (const auto &[e, c] : p.t) {
    double m = 1.0;
    for (int k = 0; k < p.nvars; ++k) for (int j = 0; j < e[static_cast<std::size_t>(k)]; ++j) m *= x[static_cast<std::size_t>(k)];
    s = s + Cx{c.re * m, c.im * m};
  }
  return s;
}

// ---- 4x4 spinor matrix of polynomials --------------------------------------------------------
struct Mat4 {
  int nvars;
  std::array<std::array<MPoly, 4>, 4> a;
  explicit Mat4(int nv) : nvars(nv) { for (auto &row : a) for (auto &e : row) e = MPoly(nv); }
};
static Mat4 matmul(const Mat4 &A, const Mat4 &B) {
  Mat4 C(A.nvars);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      MPoly s(A.nvars);
      for (int k = 0; k < 4; ++k) s = s + A.a[i][k] * B.a[k][j];
      C.a[i][j] = s;
    }
  return C;
}
static MPoly mtrace(const Mat4 &A) {
  MPoly s(A.nvars);
  for (int i = 0; i < 4; ++i) s = s + A.a[i][i];
  return s;
}
// slash p̸ = Σ_mu p^mu γ^mu, where momentum m has components = vars [4*m .. 4*m+3] (ABSTRACT basis).
static Mat4 slash(int nvars, int m) {
  Mat4 S(nvars);
  for (int mu = 0; mu < 4; ++mu) {
    MPoly comp = MPoly::var(nvars, 4 * m + mu); // p_m^mu
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const Cx g = kGamma[mu][i][j];
        if (g.re == 0 && g.im == 0) continue;
        S.a[i][j] = S.a[i][j] + MPoly::constant(nvars, g) * comp;
      }
  }
  return S;
}
// slash from explicit per-component polynomials (for momenta whose components are symbolic, e.g. the
// loop momentum whose direction carries the symbolic angles cos1/sin1/...).
static Mat4 slashC(int nvars, const std::array<MPoly, 4> &comp) {
  Mat4 S(nvars);
  for (int mu = 0; mu < 4; ++mu) {
    if (comp[static_cast<std::size_t>(mu)].t.empty()) continue;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const Cx g = kGamma[mu][i][j];
        if (g.re == 0 && g.im == 0) continue;
        S.a[i][j] = S.a[i][j] + MPoly::constant(nvars, g) * comp[static_cast<std::size_t>(mu)];
      }
  }
  return S;
}
// FRAME basis: momentum m = (symbolic magnitude var `m`) × (NUMERIC unit direction `dir`). The only
// symbolic var per momentum is its magnitude → the trace collapses to a polynomial in the magnitudes
// (≈ the sp-monomial count), the regime where the matrix-product approach can beat sp-expansion.
static Mat4 slashFrame(int nvars, int m, const std::array<double, 4> &dir) {
  Mat4 S(nvars);
  MPoly mag = MPoly::var(nvars, m); // |p_m|
  for (int mu = 0; mu < 4; ++mu) {
    if (dir[static_cast<std::size_t>(mu)] == 0.0) continue;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const Cx g = kGamma[mu][i][j];
        if (g.re == 0 && g.im == 0) continue;
        S.a[i][j] = S.a[i][j] + MPoly::constant(nvars, Cx{g.re * dir[static_cast<std::size_t>(mu)], g.im * dir[static_cast<std::size_t>(mu)]}) * mag;
      }
  }
  return S;
}

// ---- ground truth: direct NUMERIC 4x4 slash chain trace --------------------------------------
struct NMat { Cx a[4][4]; };
static NMat numSlash(const std::vector<double> &x, int m) {
  NMat S{};
  for (int mu = 0; mu < 4; ++mu) {
    const double pm = x[static_cast<std::size_t>(4 * m + mu)];
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) { const Cx g = kGamma[mu][i][j]; S.a[i][j] = S.a[i][j] + Cx{g.re * pm, g.im * pm}; }
  }
  return S;
}
static NMat nmul(const NMat &A, const NMat &B) {
  NMat C{};
  for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) { Cx s{0, 0}; for (int k = 0; k < 4; ++k) s = s + A.a[i][k] * B.a[k][j]; C.a[i][j] = s; }
  return C;
}

int main() {
  std::printf("== matrix-product prototype: polynomial-entry gamma trace ==\n");
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  int fails = 0;
  for (int K : {2, 4, 6, 8}) {
    const int nvars = 4 * K; // abstract momentum components p_m^mu
    // symbolic chain trace
    Mat4 chain = slash(nvars, 0);
    for (int m = 1; m < K; ++m) chain = matmul(chain, slash(nvars, m));
    MPoly tr = mtrace(chain);
    // numeric check at a random point
    std::vector<double> x(static_cast<std::size_t>(nvars));
    for (double &v : x) v = U(rng);
    NMat nchain = numSlash(x, 0);
    for (int m = 1; m < K; ++m) nchain = nmul(nchain, numSlash(x, m));
    Cx ntr{0, 0}; for (int i = 0; i < 4; ++i) ntr = ntr + nchain.a[i][i];
    Cx str = eval(tr, x);
    const double err = std::abs(str.re - ntr.re) + std::abs(str.im - ntr.im);
    std::printf("  K=%d  trace monomials=%zu  |sym-num|=%.2e  %s\n", K, tr.size(), err,
                err < 1e-9 ? "ok" : "FAIL");
    if (err >= 1e-9) ++fails;
  }
  // ---- FRAME basis: momenta = magnitude_var × numeric direction (the compactness regime) --------
  std::printf("\n== frame basis (momenta = |p_m| x numeric direction): trace size in magnitudes ==\n");
  for (int K : {2, 4, 6, 8, 10, 12}) {
    const int nvars = K; // one magnitude variable per momentum
    std::vector<std::array<double, 4>> dir(static_cast<std::size_t>(K));
    for (int m = 0; m < K; ++m) { // fixed pseudo-random unit directions
      std::array<double, 4> d{}; double n2 = 0;
      for (int mu = 0; mu < 4; ++mu) { d[static_cast<std::size_t>(mu)] = U(rng); n2 += d[static_cast<std::size_t>(mu)] * d[static_cast<std::size_t>(mu)]; }
      const double inv = 1.0 / std::sqrt(n2);
      for (double &v : d) v *= inv;
      dir[static_cast<std::size_t>(m)] = d;
    }
    Mat4 chain = slashFrame(nvars, 0, dir[0]);
    for (int m = 1; m < K; ++m) chain = matmul(chain, slashFrame(nvars, m, dir[static_cast<std::size_t>(m)]));
    MPoly tr = mtrace(chain);
    // expected sp-monomial count for a 2n-gamma trace = (2n-1)!! pairings (an upper bound on distinct)
    long long dfac = 1; for (int n = K - 1; n > 0; n -= 2) dfac *= n;
    std::printf("  K=%2d  frame-trace monomials=%-5zu   (sp pairings (2n-1)!! = %lld)\n", K, tr.size(), K % 2 ? 0 : dfac);
  }
  // ---- realistic: loop momentum l (symbolic ANGLES) repeating, interleaved with distinct externals.
  //      vars: 0=l1, 1=cos1, 2=sin1, 3=cos2, 4=sin2, 5=p.  l = l1*(cos1, sin1*cos2, sin1*sin2, 0).
  //      tr( l̸ p̸_1 l̸ p̸_2 ... ) with l appearing nL times — the regime of a real loop integrand. --
  std::printf("\n== realistic frame: loop momentum l (symbolic angles) repeating + externals ==\n");
  {
    const int NV = 6; // l1,cos1,sin1,cos2,sin2,p
    auto L1 = MPoly::var(NV, 0), C1 = MPoly::var(NV, 1), S1 = MPoly::var(NV, 2),
         C2 = MPoly::var(NV, 3), S2 = MPoly::var(NV, 4), P = MPoly::var(NV, 5);
    std::array<MPoly, 4> lc = {L1 * C1, L1 * S1 * C2, L1 * S1 * S2, MPoly(NV)}; // loop momentum comps
    // a few external directions (numeric), scaled by |p|
    auto extSlash = [&](double d0, double d1, double d2, double d3) {
      std::array<MPoly, 4> c = {P * MPoly::constant(NV, Cx{d0, 0}), P * MPoly::constant(NV, Cx{d1, 0}),
                                P * MPoly::constant(NV, Cx{d2, 0}), P * MPoly::constant(NV, Cx{d3, 0})};
      return slashC(NV, c);
    };
    Mat4 lS = slashC(NV, lc);
    std::array<Mat4, 3> ext = {extSlash(1, 0, 0, 0), extSlash(-0.5, 0.866, 0, 0), extSlash(-0.5, -0.866, 0, 0)};
    for (int nL : {2, 3, 4, 6}) { // even 2*nL-gamma trace: (l e)^nL, l repeats nL times
      Mat4 chain = lS;
      for (int r = 0; r < nL; ++r) {
        chain = matmul(chain, ext[static_cast<std::size_t>(r % 3)]);
        if (r < nL - 1) chain = matmul(chain, lS);
      }
      MPoly tr = mtrace(chain);
      std::printf("  l repeats %d times (%d-gamma trace):  frame monomials in {l1,cos1,sin1,cos2,sin2,p} = %zu\n",
                  nL, 2 * nL, tr.size());
    }
    std::printf("  (vs reduce_product: the same structures are where it emits ~10^9 sp-monomials)\n");
  }
  std::printf(fails == 0 ? "\nALL TESTS PASSED\n" : "\nTESTS FAILED\n");
  return fails ? 1 : 0;
}
