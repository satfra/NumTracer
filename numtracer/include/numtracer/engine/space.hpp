/// @file space.hpp
/// @brief @ref numtracer::Frame — the kinematic frame a network is contracted on: the scalar symbols,
///        the momenta's four components in those symbols, and the projector denominators `k²`.
///
/// Everything a contraction needs besides the network itself lives here, so a contraction is one
/// call, `F.trace(chain, net)`, and an evaluation is `F.eval(poly, F.at(...))`:
///
/// ```cpp
/// nt::Frame F;
/// auto P = F.symbol("p"), L = F.symbol("l");
/// auto [C, S] = F.angle("theta");                // S = sqrt(1 - C^2)
/// auto p = F.momentum(P, 0, 0, 0);
/// auto l = F.momentum(L * C, L * S, 0, 0);
/// auto [mu, nu] = F.indices<2>();
/// nt::Poly T = F.trace({nt::slash(p), nt::gamma(mu), nt::slash(l - p), nt::gamma(nu)}, nt::projT(mu, nu, l));
/// double v = F.eval(T, F.at(1.3, 0.8, 0.6)).re;
/// ```
///
/// A frame is built in two stages. First declare its symbols (`symbol`, `angle`); the first time a
/// polynomial is made from them — a momentum, an arithmetic expression — the symbol list is FROZEN,
/// because every polynomial of one frame must share it. Declaring a symbol after that throws.
///
/// Generated code uses the same object: it constructs the frame with its symbol list up front, fills
/// the component table with @ref Frame::set_component, registers the denominators of all its nets
/// once with @ref Frame::add_denominators, and then contracts in parallel through the `const`
/// members.
#pragma once

#include "numtracer/core/config.hpp" // NT_THROW

#include <stdexcept>
#include <string>

namespace numtracer::inline engine
{

  struct Space {
    int dim = 4;
    int spinor = 4;
  };

  inline void validate(const Space &space)
  {
    if (space.dim < 1 || space.dim > 4)
      NT_THROW(std::invalid_argument, ("Space: invalid dimension: " + std::to_string(space.dim)).c_str());
    if ((space.spinor != 2 && space.spinor != 4) || (space.dim == 4 && space.spinor != 4))
      NT_THROW(std::invalid_argument, ("Space: invalid spinor dimension: " + std::to_string(space.spinor)).c_str());
  }

} // namespace numtracer::inline engine