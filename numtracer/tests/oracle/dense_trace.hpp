/// @file dense_trace.hpp
/// @brief Test-only chiral Dirac gamma-trace oracle (2×2 block reduction).
///
/// An independent, brute-structure evaluator for the all-slash Dirac chain, kept in the
/// test tree (`tests/oracle/`) as ground truth for the numeric contraction engine's unit
/// gate (`test_contract`). It is not part of the shipped engine.
///
/// A chain
/// `tr(p̸₁ … p̸_K)` evaluated by multiplying slashed-momentum matrices and taking
/// the matrix trace. Instead of a brute-force 4×4 matmul, it
/// exploits the **chiral block structure** of the typed-out Weyl basis: every
/// slashed momentum is block-antidiagonal,
/// @f[ p̸ = \begin{pmatrix} 0 & P \\ Q & 0 \end{pmatrix}, @f]
/// with 2×2 complex blocks `P` (upper-right) and `Q` (lower-left). A product of
/// such matrices stays "block-monomial": after an even number of factors it is
/// block-diagonal `[[U,0],[0,V]]`, after an odd number block-antidiagonal again.
/// So the whole 4×4 chain reduces to a fold of **2×2** products — a quarter of the
/// arithmetic, each 2×2 multiply a tiny straight-line kernel the compiler
/// vectorises (measured ~3–6× faster than a brute 4×4 chain). The block
/// structure also makes `tr(odd chain) = 0` fall out for free: an odd chain stays
/// block-antidiagonal, so its diagonal — hence its trace — is zero.
///
/// (A general structurally-sparse 4×4 matmul was tried and is *slower* than the
/// brute path: at 4×4 the matmul is vectorisation-bound, and zero-skipping — even
/// branch-free at compile time — serialises the work. The block reduction wins
/// because it does genuinely less work in a still-vectorisable shape.)
#pragma once

#include "numtracer/core/matrix.hpp"      // Matrix
#include "numtracer/dirac/dirac_data.hpp" // kGamma

#include <complex>

namespace numtracer::dirac {

/// @brief Spinor dimension (4 in 4D).
inline constexpr int kSpinDim = 4;

namespace dirac_detail {
/// @brief 2×2 complex matrix product (fully unrolled).
/// @param x Left matrix.
/// @param y Right matrix.
/// @return `x · y`.
inline Matrix<2> mul2(const Matrix<2> &x, const Matrix<2> &y) {
  return Matrix<2>{{x(0, 0) * y(0, 0) + x(0, 1) * y(1, 0), x(0, 0) * y(0, 1) + x(0, 1) * y(1, 1),
                    x(1, 0) * y(0, 0) + x(1, 1) * y(1, 0), x(1, 0) * y(0, 1) + x(1, 1) * y(1, 1)}};
}
/// @brief 2×2 trace `a00 + a11`.
/// @param x The matrix.
/// @return The trace.
inline std::complex<double> tr2(const Matrix<2> &x) { return x(0, 0) + x(1, 1); }

/// @brief The two nonzero 2×2 blocks of a dressed slashed momentum `z·p̸`.
///
/// `P` is the upper-right block (`kGamma[μ][r][c+2]`), `Q` the lower-left
/// (`kGamma[μ][r+2][c]`), summed over the momentum components.
/// @param p The 4-momentum.
/// @param z The real dressing factor.
/// @param P Out: the upper-right block.
/// @param Q Out: the lower-left block.
inline void slashed_blocks(const double *p, double z, Matrix<2> &P, Matrix<2> &Q) {
  for (int mu = 0; mu < kSpinDim; ++mu) {
    const double pz = p[mu] * z;
    for (int r = 0; r < 2; ++r)
      for (int c = 0; c < 2; ++c) {
        const Cx gamma_upper = kGamma[mu][r][c + 2];
        P(r, c) = P(r, c) + std::complex<double>{pz * gamma_upper.re, pz * gamma_upper.im};
        const Cx gamma_lower = kGamma[mu][r + 2][c];
        Q(r, c) = Q(r, c) + std::complex<double>{pz * gamma_lower.re, pz * gamma_lower.im};
      }
  }
}
} // namespace dirac_detail

/// @brief Dense Dirac trace `tr(z₁ p̸₁ · … · z_K p̸_K)` via the chiral 2×2 blocks.
///
/// Folds the chain in block form: the running product is two 2×2 blocks whose
/// roles (diagonal vs antidiagonal) alternate with each factor. Each leg costs two
/// 2×2 products; an odd chain returns 0 structurally.
/// @tparam K The chain length (compile-time).
/// @param P The K momenta (`P[i]` is the i-th 4-momentum).
/// @param z The per-leg real dressings (length K), or `nullptr` for all 1.
/// @return The complex trace.
template <int K> std::complex<double> chiral_gamma_trace(const double (*P)[4], const double *z = nullptr) {
  Matrix<2> m0, m1; // antidiagonal blocks of leg 0: m0 = P₀ (upper), m1 = Q₀ (lower)
  dirac_detail::slashed_blocks(P[0], z ? z[0] : 1.0, m0, m1);
  bool antidiag = true; // after 1 factor the product is block-antidiagonal
  for (int i = 1; i < K; ++i) {
    Matrix<2> Pi, Qi;
    dirac_detail::slashed_blocks(P[i], z ? z[i] : 1.0, Pi, Qi);
    if (antidiag) { // [[0,m0],[m1,0]] · [[0,Pi],[Qi,0]] = [[m0·Qi,0],[0,m1·Pi]]
      m0 = dirac_detail::mul2(m0, Qi);
      m1 = dirac_detail::mul2(m1, Pi);
    } else { // [[m0,0],[0,m1]] · [[0,Pi],[Qi,0]] = [[0,m0·Pi],[m1·Qi,0]]
      m0 = dirac_detail::mul2(m0, Pi);
      m1 = dirac_detail::mul2(m1, Qi);
    }
    antidiag = !antidiag;
  }
  if (antidiag) return std::complex<double>{0, 0}; // odd chain: zero diagonal
  return dirac_detail::tr2(m0) + dirac_detail::tr2(m1);
}

} // namespace numtracer::dirac
