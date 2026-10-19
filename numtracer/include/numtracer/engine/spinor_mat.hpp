/// @file spinor_mat.hpp
/// @brief A 4×4 spinor matrix whose entries are @ref numtracer::poly::Poly, plus the γ/slash
///        builders the numeric Dirac trace multiplies.
///
/// The numeric contraction backend contracts a Dirac chain as a product of 4×4 matrices whose
/// entries are polynomials in the frame's scalar symbols (@ref numtracer::poly::Poly): the γ
/// matrices stay numeric (constant entries) while only the user's momentum data enters as
/// polynomial coefficients. This header holds the leaf builders of that `Matrix<4, Poly>`:
/// `gammaC` (a numeric γ^μ), `cmatC` (the charge-conjugation matrix C) and `slashC` (a slashed
/// momentum from explicit per-component polynomials). Product and trace are the generic
/// `matmul`/`mtrace` of `core/cmat.hpp`.
#pragma once

#include "numtracer/core/cmat.hpp"        // Matrix
#include "numtracer/dirac/dirac_data.hpp" // kGamma, kC (Euclidean Weyl)
#include "numtracer/poly/mpoly.hpp"       // Poly

#include <array>
#include <utility> // std::move

namespace numtracer::inline engine
{

  /// @brief An N×N matrix of zero polynomials over `nsym` symbols.
  ///
  /// `Matrix<N, Poly>{}` holds *default* polynomials (`nsym == 0`); this gives every entry the
  /// frame's symbol count instead, which the γ/C/slash builders below start from.
  template <int N> Matrix<N, Poly> zeroPolyMatrix(int nsym)
  {
    Matrix<N, Poly> m;
    for (auto &e : m.data)
      e = PolyFactory::zero(nsym);
    return m;
  }

  /// @brief The numeric γ matrix for a concrete Lorentz index `mu` (a constant 4×4 of `Poly`).
  inline Matrix<4, Poly> gammaC(int nsym, int mu)
  {
    Matrix<4, Poly> S = zeroPolyMatrix<4>(nsym);
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const Cx g = numtracer::dirac::kGamma[mu][i][j];
        if (g.re == 0 && g.im == 0) continue;
        S(i, j) = PolyFactory::constant(nsym, g);
      }
    return S;
  }

  /// @brief The numeric charge-conjugation matrix `C` (a constant 4×4 of `Poly`).
  ///
  /// Reads @ref numtracer::dirac::kC, which is the single source of truth for `C`; the fold in
  /// `trace_dirac_loop` takes this matrix's DIAGONAL Weyl blocks (C is block-diagonal, unlike γ^μ).
  inline Matrix<4, Poly> cmatC(int nsym)
  {
    Matrix<4, Poly> S = zeroPolyMatrix<4>(nsym);
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const Cx c = numtracer::dirac::kC[i][j];
        if (c.re == 0 && c.im == 0) continue;
        S(i, j) = PolyFactory::constant(nsym, c);
      }
    return S;
  }

  /// @brief Slashed momentum `p̸ = Σ_μ comp[μ] γ^μ` from explicit per-component polynomials.
  inline Matrix<4, Poly> slashC(int nsym, const std::array<Poly, 4> &comp)
  {
    Matrix<4, Poly> S = zeroPolyMatrix<4>(nsym);
    for (int mu = 0; mu < 4; ++mu) {
      if (comp[mu].terms.empty()) continue;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          const Cx g = numtracer::dirac::kGamma[mu][i][j];
          if (g.re == 0 && g.im == 0) continue;
          // `scaled` instead of `constant(g) * comp[mu]`: bit-identical (same operand and monomial
          // order) without the scratch sort — this runs per Slash token inside the Dirac fold.
          S(i, j) = std::move(S(i, j)) + PolyFactory::scaled(nsym, comp[mu], g);
        }
    }
    return S;
  }

} // namespace numtracer::inline engine
