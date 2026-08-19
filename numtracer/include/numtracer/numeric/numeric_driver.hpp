/// @file numeric_driver.hpp
/// @brief Render an @ref MPoly (a projector denominator `k²`, or any component expression) as a C++
///        expression in the user symbol names, for the generator's `FillFormulas` — the
///        `inv(atom) = 1/k²` and `var(k) = <name>` slots.
///
/// This file used to also hold a `run_numeric` driver that contracted and lowered every net in one
/// serial loop. Codegen.m emits the parallel `contract_traces` + `fold_groups_streaming` pair
/// (numeric/trace_fold.hpp) instead, and had done for long enough that nothing called the driver;
/// it and its `NumericNet` net type are gone.
#pragma once

#include "numtracer/numeric/numeric_contract.hpp"
#include "numtracer/core/config.hpp" // NT_THROW (exception-optional guard for -fno-exceptions builds)

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace numtracer::numeric {

/// @brief Render an @ref MPoly as a C++ expression in @p symNames (real coefficients only — used
///        for projector denominators `k²` and component expressions, which carry no imaginary part
///        and no inverse atoms). Powers are emitted as repeated multiplication.
inline std::string mpoly_to_cpp(const MPoly &p, const std::vector<std::string> &symNames) {
  if (p.terms.empty()) return "0.0";
  if ((int)symNames.size() < p.nsym)
    NT_THROW(std::runtime_error, "mpoly_to_cpp: symNames shorter than the polynomial's symbol count");
  std::ostringstream os;
  os.setf(std::ios::scientific);
  os.precision(17);
  // an imaginary part above this is a genuine complex coefficient (not round-off): this renderer is
  // for real-only expressions (projector denominators / component formulas), so it rejects those.
  constexpr double kRealCoeffTol = 1e-12;
  bool first = true;
  for (const auto &[m, c] : p.terms) {
    if (std::abs(c.im) > kRealCoeffTol)
      NT_THROW(std::runtime_error, "mpoly_to_cpp: complex coefficient where a real expression was expected");
    if (!m.atoms.empty())
      NT_THROW(std::runtime_error, "mpoly_to_cpp: monomial carries an inverse atom (not a plain expression)");
    const double v = c.re;
    if (!first) os << (v < 0 ? " - " : " + ");
    else if (v < 0) os << "-";
    const double av = v < 0 ? -v : v;
    os << "(" << av << ")";
    for (int k = 0; k < p.nsym; ++k)
      for (int e = 0; e < m.e[k]; ++e) os << "*" << symNames[k];
    first = false;
  }
  return os.str();
}

} // namespace numtracer::numeric
