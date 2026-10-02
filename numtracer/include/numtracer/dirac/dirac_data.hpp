/// @file dirac_data.hpp
/// @brief Euclidean Dirac gamma matrices as typed-out `constexpr` tables (Weyl basis).
///
/// 4D Euclidean, Hermitian gammas in the chiral/Weyl representation:
/// @f[
///   \gamma^k = \begin{pmatrix} 0 & -i\sigma^k \\ i\sigma^k & 0 \end{pmatrix}
///   \quad (k = 1,2,3), \qquad
///   \gamma^4 = \begin{pmatrix} 0 & I \\ I & 0 \end{pmatrix}.
/// @f]
/// They satisfy the Euclidean Clifford algebra
/// @f$\{\gamma^\mu,\gamma^\nu\} = 2\,\delta^{\mu\nu} I@f$, and
/// @f$\gamma_5 = \mathrm{diag}(1,1,-1,-1)@f$ anticommutes with all four.
///
/// Index conventions: `mu = 0..3` corresponds to @f$\gamma^1..\gamma^4@f$; spinor
/// indices `i, j = 0..3`. The values are stored as plain `constexpr` tables, read by
/// `numeric/spinor_mat.hpp` (and the test-side dense oracle, `tests/oracle/dtensor.hpp`), and
/// are cross-checked against the Clifford algebra in `tests/test_numeric_contract.cpp`
/// (sections A/F for the gammas, H for `kC` below) -- that test is the safety net for
/// these tables, which nothing else would notice were wrong.
#pragma once

#include "numtracer/core/cx.hpp" // Cx

namespace numtracer::dirac {

/// @brief The gamma matrices: `kGamma[mu][i][j]` is @f$(\gamma^\mu)_{ij}@f$.
///
/// The imaginary unit is `Cx{0,1}` and `-i` is `Cx{0,-1}`. `mu` runs `0..3` over
/// @f$\gamma^1..\gamma^4@f$; `i, j` are spinor indices `0..3`.
inline constexpr Cx kGamma[4][4][4] = {
    // gamma^1 = [[0,-i sigma1],[i sigma1,0]]
    {{Cx{0, 0}, Cx{0, 0}, Cx{0, 0}, Cx{0, -1}},
     {Cx{0, 0}, Cx{0, 0}, Cx{0, -1}, Cx{0, 0}},
     {Cx{0, 0}, Cx{0, 1}, Cx{0, 0}, Cx{0, 0}},
     {Cx{0, 1}, Cx{0, 0}, Cx{0, 0}, Cx{0, 0}}},
    // gamma^2 = [[0,-i sigma2],[i sigma2,0]]
    {{Cx{0, 0}, Cx{0, 0}, Cx{0, 0}, Cx{-1, 0}},
     {Cx{0, 0}, Cx{0, 0}, Cx{1, 0}, Cx{0, 0}},
     {Cx{0, 0}, Cx{1, 0}, Cx{0, 0}, Cx{0, 0}},
     {Cx{-1, 0}, Cx{0, 0}, Cx{0, 0}, Cx{0, 0}}},
    // gamma^3 = [[0,-i sigma3],[i sigma3,0]]
    {{Cx{0, 0}, Cx{0, 0}, Cx{0, -1}, Cx{0, 0}},
     {Cx{0, 0}, Cx{0, 0}, Cx{0, 0}, Cx{0, 1}},
     {Cx{0, 1}, Cx{0, 0}, Cx{0, 0}, Cx{0, 0}},
     {Cx{0, 0}, Cx{0, -1}, Cx{0, 0}, Cx{0, 0}}},
    // gamma^4 = [[0,I],[I,0]]
    {{Cx{0, 0}, Cx{0, 0}, Cx{1, 0}, Cx{0, 0}},
     {Cx{0, 0}, Cx{0, 0}, Cx{0, 0}, Cx{1, 0}},
     {Cx{1, 0}, Cx{0, 0}, Cx{0, 0}, Cx{0, 0}},
     {Cx{0, 0}, Cx{1, 0}, Cx{0, 0}, Cx{0, 0}}}};

/// @brief @f$\gamma_5 = \mathrm{diag}(1,1,-1,-1)@f$ in the chiral basis.
inline constexpr Cx kGamma5[4][4] = {{Cx{1, 0}, Cx{0, 0}, Cx{0, 0}, Cx{0, 0}},
                                     {Cx{0, 0}, Cx{1, 0}, Cx{0, 0}, Cx{0, 0}},
                                     {Cx{0, 0}, Cx{0, 0}, Cx{-1, 0}, Cx{0, 0}},
                                     {Cx{0, 0}, Cx{0, 0}, Cx{0, 0}, Cx{-1, 0}}};

/// @brief The charge-conjugation matrix @f$C = \gamma^2\gamma^4@f$ in this (Euclidean, Weyl) basis.
///
/// Defined so that it satisfies the standard identities, all verified against these tables in
/// `tests/test_numeric_contract.cpp` section H:
/// @f[
///   C^T = C^{-1} = -C^\dagger = -C, \qquad
///   C\gamma_\mu^T C^{-1} = -\gamma_\mu, \qquad
///   C\gamma_5^T C^{-1} = +\gamma_5 .
/// @f]
/// This matches the convention of arXiv:2606.23772 Eq. (46), so a diquark vertex written there
/// transcribes here unchanged.
///
/// Two properties matter to the engine. (1) `C` is block-DIAGONAL in the Weyl split — it is a
/// product of two block-antidiagonal factors — so, exactly like @f$\gamma_5@f$, it does not flip
/// the antidiagonal trace parity and needs no full 2x2 multiply; its two blocks are signed
/// permutations. (2) @f$C\gamma_5@f$ is ANTISYMMETRIC, which is what makes the scalar diquark
/// bilinear @f$q^T C\gamma_5 q@f$ non-vanishing once the antisymmetric colour and flavour factors
/// are included.
inline constexpr Cx kC[4][4] = {{Cx{0, 0}, Cx{-1, 0}, Cx{0, 0}, Cx{0, 0}},
                                {Cx{1, 0}, Cx{0, 0}, Cx{0, 0}, Cx{0, 0}},
                                {Cx{0, 0}, Cx{0, 0}, Cx{0, 0}, Cx{1, 0}},
                                {Cx{0, 0}, Cx{0, 0}, Cx{-1, 0}, Cx{0, 0}}};

/// @brief Look up a gamma-matrix entry.
/// @param mu The matrix index `0..3` (selecting @f$\gamma^1..\gamma^4@f$).
/// @param i The row (spinor) index `0..3`.
/// @param j The column (spinor) index `0..3`.
/// @return @f$(\gamma^\mu)_{ij}@f$.
constexpr Cx gamma_entry(int mu, int i, int j) { return kGamma[mu][i][j]; }
/// @brief Look up a @f$\gamma_5@f$ entry.
/// @param i The row (spinor) index `0..3`.
/// @param j The column (spinor) index `0..3`.
/// @return @f$(\gamma_5)_{ij}@f$.
constexpr Cx gamma5_entry(int i, int j) { return kGamma5[i][j]; }
/// @brief Look up a charge-conjugation-matrix entry.
/// @param i The row (spinor) index `0..3`.
/// @param j The column (spinor) index `0..3`.
/// @return @f$C_{ij}@f$.
constexpr Cx c_entry(int i, int j) { return kC[i][j]; }

} // namespace numtracer::dirac
