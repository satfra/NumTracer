/// @file numeric_driver.hpp
/// @brief Build-time driver for the numeric (matrix-product) backend: contract every diagram with
///        @ref numeric_value and lower it into ONE shared fundamental-symbol env, returning a
///        per-net @ref network::GenProg plus the env (so the generator emits `fill` / `trN` / the kernel
///        assembly exactly as the inv path does). No `reduce`/`rebase`/`ibp` — the contraction is
///        already a small polynomial — so this is serial and negligible in RAM/time.
///
/// Also provides @ref mpoly_to_cpp: render an @ref MPoly (a projector denominator `k²`, or any
/// component expression) as a C++ expression in the user symbol names, for the generator's
/// `FillFormulas` (the `inv(atom) = 1/k²` and `var(k) = <name>` slots).
#pragma once

#include "numtracer/numeric/numeric_contract.hpp"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace numtracer::numeric {

/// @brief One diagram in the numeric backend: a Dirac chain (possibly empty) and a Lorentz network.
struct NumericNet {
  network::DiracNet dirac;
  NNet lorentz;
};

/// @brief Contract and lower every net into the shared env @p g; returns one program per net.
inline std::vector<network::GenProg> run_numeric(int nsym, const std::vector<NumericNet> &nets,
                                             const std::vector<std::array<MPoly, 4>> &comp,
                                             const std::vector<MPoly> &atomDen, network::GlobalEnv &g) {
  std::vector<network::GenProg> progs;
  progs.reserve(nets.size());
  for (const NumericNet &n : nets) {
    MPoly mp = numeric_value(nsym, n.dirac, n.lorentz, comp, atomDen);
    progs.push_back(to_genprog(mp, g));
  }
  return progs;
}

/// @brief Render an @ref MPoly as a C++ expression in @p symNames (real coefficients only — used
///        for projector denominators `k²` and component expressions, which carry no imaginary part
///        and no inverse atoms). Powers are emitted as repeated multiplication.
inline std::string mpoly_to_cpp(const MPoly &p, const std::vector<std::string> &symNames) {
  if (p.t.empty()) return "0.0";
  std::ostringstream os;
  os.setf(std::ios::scientific);
  os.precision(17);
  bool first = true;
  for (const auto &[m, c] : p.t) {
    if (std::abs(c.im) > 1e-12)
      throw std::runtime_error("mpoly_to_cpp: complex coefficient where a real expression was expected");
    if (!m.atoms.empty())
      throw std::runtime_error("mpoly_to_cpp: monomial carries an inverse atom (not a plain expression)");
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
