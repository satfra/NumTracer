/// @file numeric_driver.hpp
/// @brief Render an @ref Poly (a projector denominator `k²`, or any component expression) as a C++
///        expression in the user symbol names, for the generator's `FillFormulas` — the
///        `inv(atom) = 1/k²` and `var(k) = <name>` slots, including the (Re, Im) halves of a COMPLEX
///        `1/k²` (a finite-density frame; see @ref inv_fill_cpp).
///
/// The contraction driver itself is `contract_traces` + `fold_groups_streaming`
/// (numeric/trace_fold.hpp), emitted by mathematica/CodegenGenerator.m.
#pragma once

#include "numtracer/numeric/mpoly.hpp"
#include "numtracer/core/config.hpp" // NT_THROW (exception-optional guard for -fno-exceptions builds)
#include "numtracer/codegen/precision.hpp"
#include "numtracer/codegen/gen.hpp" // network::GlobalEnv (mark_complex_atoms)

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace numtracer::inline numeric {

/// @brief An imaginary part above this is a genuine complex coefficient (not round-off).
inline constexpr double kRealCoeffTol = 1e-12;

/// @brief Which half of a (possibly complex) polynomial @ref poly_render emits.
enum class PolyPart {
  RealStrict, ///< the polynomial must BE real; an imaginary coefficient throws (@ref poly_to_cpp)
  Re,         ///< the real half
  Im          ///< the imaginary half (itself a real-coefficient polynomial)
};

/// @brief Render one half of an @ref Poly as a real C++ expression in @p symNames. Powers are
///        emitted as repeated multiplication; an inverse atom is rejected.
inline std::string poly_render(const Poly &p, const std::vector<std::string> &symNames, PolyPart part) {
  const char *zero = codegen::emit_single() ? "0.f" : "0.0";
  if (p.terms.empty()) return zero;
  if ((int)symNames.size() < p.nsym)
    NT_THROW(std::runtime_error, "poly_to_cpp: symNames shorter than the polynomial's symbol count");
  std::ostringstream os;
  os.setf(std::ios::scientific);
  os.precision(17);
  bool first = true;
  for (const auto &[m, c] : p.terms) {
    if (part == PolyPart::RealStrict && std::abs(c.im) > kRealCoeffTol)
      NT_THROW(std::runtime_error, "poly_to_cpp: complex coefficient where a real expression was expected");
    if (!m.atoms.empty())
      NT_THROW(std::runtime_error, "poly_to_cpp: monomial carries an inverse atom (not a plain expression)");
    const double v = part == PolyPart::Im ? c.im : c.re;
    // a half drops the terms living entirely in the other half; RealStrict keeps its exact output
    if (part != PolyPart::RealStrict && v == 0.0) continue;
    if (!first) os << (v < 0 ? " - " : " + ");
    else if (v < 0) os << "-";
    const double av = v < 0 ? -v : v;
    os << "(";
    if (codegen::emit_single())
      os << codegen::float_literal(av);
    else
      os << av;
    os << ")";
    for (int k = 0; k < p.nsym; ++k)
      for (int e = 0; e < m.e[k]; ++e) os << "*" << symNames[k];
    first = false;
  }
  if (first) return zero;
  return os.str();
}

/// @brief Render an @ref Poly as a C++ expression in @p symNames (real coefficients only — used
///        for projector denominators `k²` and component expressions, which carry no imaginary part
///        and no inverse atoms). Powers are emitted as repeated multiplication.
inline std::string poly_to_cpp(const Poly &p, const std::vector<std::string> &symNames) {
  return poly_render(p, symNames, PolyPart::RealStrict);
}

/// @brief Does @p p carry a genuine imaginary part?
inline bool poly_complex(const Poly &p) {
  for (const auto &[m, c] : p.terms)
    if (std::abs(c.im) > kRealCoeffTol) return true;
  return false;
}

/// @brief The `fill` formula of an inverse atom `1/k²`, or of its imaginary half.
///
/// A complex `k² = A + i B` (an internal line mixing the loop with a quark leg of frequency
/// p0 - i mu) splits into the real slot pair `Re = A/(A²+B²)`, `Im = -B/(A²+B²)`. A real `k²`
/// renders as `1.0/(k²)` (`1.f/(k²)` single), the spelling the generator always used.
inline std::string inv_fill_cpp(const Poly &p, const std::vector<std::string> &symNames, bool imagPart) {
  const char *one = codegen::emit_single() ? "1.f" : "1.0";
  if (!poly_complex(p))
    return imagPart ? std::string(codegen::emit_single() ? "0.f" : "0.0") : one + ("/(" + poly_to_cpp(p, symNames) + ")");
  const std::string A = poly_render(p, symNames, PolyPart::Re);
  const std::string B = poly_render(p, symNames, PolyPart::Im);
  const std::string den = "((" + A + ")*(" + A + ") + (" + B + ")*(" + B + "))";
  return imagPart ? ("-(" + B + ")/" + den) : ("(" + A + ")/" + den);
}

/// @brief Flag the inverse atoms whose denominator is complex, so the lowering splits them into
///        `(Re, Im)` slot pairs. Call once, after the denominators are final and BEFORE any trace is
///        lowered: the flags decide how every monomial carrying the atom is interned. Leaves the env
///        untouched when every denominator is real.
inline void mark_complex_atoms(network::GlobalEnv &g, const std::vector<Poly> &atomDen) {
  for (std::size_t a = 0; a < atomDen.size(); ++a)
    if (poly_complex(atomDen[a])) {
      if (g.cplxAtom.empty()) g.cplxAtom.assign(atomDen.size(), 0);
      g.cplxAtom[a] = 1;
    }
}

} // namespace numtracer::numeric
