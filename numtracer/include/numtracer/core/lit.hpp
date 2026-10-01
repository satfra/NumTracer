/// @file lit.hpp
/// @brief Compile-time complex-constant carrier for the code generator.
///
/// `Lit<C>` is an empty type whose sole job is to carry a `constexpr` @ref numtracer::Cx value
/// `C` as a non-type template parameter. The numeric backend's generator emits scalar
/// coefficients as `numtracer::Lit<Cx{...}>` type arguments to the generated `litco`/`sc<L>`
/// helpers (see `mathematica/CodegenGenerator.m`), which recover the value via a trait specialisation.
/// It is the only surviving piece of the former compile-time expression algebra — a plain
/// constant wrapper, not part of the (runtime) lowering machinery.
#pragma once

#include "numtracer/core/cx.hpp" // Cx (constexpr complex NTTP)

namespace numtracer
{

  /// @brief A compile-time complex constant carried as a non-type template parameter.
  /// @tparam C The constant value.
  template <Cx C> struct Lit {};

} // namespace numtracer
