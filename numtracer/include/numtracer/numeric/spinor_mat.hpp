/// @file spinor_mat.hpp
/// @brief A 4×4 spinor matrix whose entries are @ref numtracer::numeric::MPoly, plus the γ/slash
///        builders the numeric Dirac trace multiplies.
///
/// The numeric contraction backend contracts a Dirac chain as a product of 4×4 matrices whose
/// entries are polynomials in the frame's scalar symbols (@ref numtracer::numeric::MPoly): the γ
/// matrices stay numeric (constant entries) while only the user's momentum data enters as
/// polynomial coefficients. This header holds that matrix type (`Mat4`), its product/trace, and
/// the two leaf builders — `gammaC` (a numeric γ^μ) and `slashC` (a slashed momentum from explicit
/// per-component polynomials). It depends on `mpoly.hpp` (for `MPoly`) and `dirac/dirac_data.hpp`
/// (for the typed-out Weyl γ table `kGamma`); `Mat4` cannot live in `core/` because its entries are
/// `MPoly`, a numeric-backend type.
#pragma once

#include "numtracer/dirac/dirac_data.hpp" // kGamma (Euclidean Weyl)
#include "numtracer/numeric/mpoly.hpp"    // MPoly

#include <array>

namespace numtracer::numeric
{

  /// @brief A 4×4 spinor matrix whose entries are @ref MPoly (numeric γ, symbolic momentum data).
  struct Mat4 {
    int nsym;
    std::array<std::array<MPoly, 4>, 4> a;
    explicit Mat4(int ns) : nsym(ns)
    {
      for (auto &row : a)
        for (auto &e : row)
          e = MPoly(ns);
    }
  };

  inline Mat4 matmul(const Mat4 &A, const Mat4 &B)
  {
    Mat4 C(A.nsym);
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        MPoly s(A.nsym);
        for (int k = 0; k < 4; ++k)
          s = s + A.a[i][k] * B.a[k][j];
        C.a[i][j] = std::move(s);
      }
    return C;
  }
  inline MPoly mtrace(const Mat4 &A)
  {
    MPoly s(A.nsym);
    for (int i = 0; i < 4; ++i)
      s = s + A.a[i][i];
    return s;
  }

  /// @brief The numeric γ matrix for a concrete Lorentz index `mu` (a constant 4×4 of `MPoly`).
  inline Mat4 gammaC(int nsym, int mu)
  {
    Mat4 S(nsym);
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        const Cx g = numtracer::dirac::kGamma[mu][i][j];
        if (g.re == 0 && g.im == 0) continue;
        S.a[i][j] = MPoly::constant(nsym, g);
      }
    return S;
  }

  /// @brief Slashed momentum `p̸ = Σ_μ comp[μ] γ^μ` from explicit per-component polynomials.
  inline Mat4 slashC(int nsym, const std::array<MPoly, 4> &comp)
  {
    Mat4 S(nsym);
    for (int mu = 0; mu < 4; ++mu) {
      if (comp[mu].t.empty()) continue;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          const Cx g = numtracer::dirac::kGamma[mu][i][j];
          if (g.re == 0 && g.im == 0) continue;
          S.a[i][j] = S.a[i][j] + MPoly::constant(nsym, g) * comp[mu];
        }
    }
    return S;
  }

} // namespace numtracer::numeric
