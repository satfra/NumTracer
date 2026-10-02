/// @file core/hash.hpp
/// @brief The 64-bit hash mixers the engine's open-addressed indices are built on.
///
/// Pulls in nothing but `<cstdint>`, so the light `codegen/real_cse.hpp` can share it without
/// depending on the heavier `network/network.hpp`.
///
/// Nothing about the emitted code depends on the choice of hash. Both indices
/// (@ref numtracer::network::GlobalEnv in `codegen/gen.hpp`, @ref numtracer::network::rdetail::RBuilder
/// in `codegen/real_cse.hpp`) return
/// *insertion* positions and compare candidates by equality, so the hash picks buckets only — any
/// mixer with the same signature would emit byte-identical kernels.
#pragma once

#include <cstdint>

namespace numtracer
{

  /// @brief splitmix64 finaliser — scrambles a 64-bit value. `constexpr`, so it is usable from the
  ///        `constexpr` builders in `codegen/real_cse.hpp`.
  constexpr std::uint64_t splitmix64_finalise(std::uint64_t x)
  {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
  }

  /// @brief boost-style combine of two 64-bit hashes (order-sensitive).
  constexpr std::uint64_t hash_combine(std::uint64_t a, std::uint64_t b)
  {
    return splitmix64_finalise(a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2)));
  }

} // namespace numtracer
