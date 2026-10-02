/// @file tolerances.hpp
/// @brief The numeric round-off tolerances of the build-time contraction, in one place.
///
/// Changing any of these changes which terms survive, hence the emitted kernels. Validate a change
/// against the INTEGRATED numeric-vs-FORM error, not a pointwise round-off floor.
#pragma once

namespace numtracer::numeric
{
  /// @brief Relative noise-prune tolerance: a monomial whose |coefficient| is below this fraction of
  ///        the largest coefficient is round-off from the numeric frame (a ~10-order gap separates it
  ///        from physics), so @ref to_genprog drops it.
  inline constexpr double kNoisePruneRelTol = 1e-9;

  /// @brief Relative tolerance for "the division remainder vanishes" in @ref divThroughPolyAtoms,
  ///        against the dividend's largest coefficient. A different test from the noise prune, but
  ///        deliberately on the same scale: an exact cancellation lands at ~1e-16 relative, and the
  ///        polynomial already carries round-off up to the prune threshold.
  inline constexpr double kPolyDivRelTol = kNoisePruneRelTol;

  /// @brief Absolute tolerance below which an SU(N) colour/flavour factor component snaps to 0.
  ///        These factors are exact rationals (× generator traces), so anything this small is the
  ///        √3 round-off the generator-table arithmetic leaves behind, well below any genuine value.
  inline constexpr double kZeroSnapTol = 1e-9;
} // namespace numtracer::numeric
