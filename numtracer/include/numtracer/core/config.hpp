/// @file core/config.hpp
/// @brief Centralised compile-time / code-generation tunables.
///
/// These constants bound how the symbolic engine lowers an expression to executable code and
/// how the build-time generator schedules its work. They are gathered here (rather than
/// scattered through the headers that use them) so the trade-offs live in one place. The
/// defaults are chosen empirically as good general settings; each is overridable, so a
/// translation unit with unusually large or small expressions can retune without editing the
/// library. None of them change *results* — only the compile-time/run-time cost profile.
#pragma once

#include <cstddef>

// ---- real-SSA lowering tiers (numtracer::expr::eval_real) ----------------------------------
// The real lowering picks one of three strategies by SSA size: fully-unrolled straight-line
// code (fastest at run time, super-linear compile cost), `noinline`-chunked straight-line code
// (a middle tier that bounds each compiled function), or a flat interpreter (constant compile
// cost, slowest per eval). The two cutoffs below select the tier; the block size sets the
// chunk granularity. All three are overridable per-TU with `-D` (e.g. `-DNT_UNROLL_MAX=`),
// which wins over these defaults because the command-line macro is seen first.

/// @brief Largest real-SSA size still emitted as fully-unrolled straight-line code.
#ifndef NT_UNROLL_MAX
#define NT_UNROLL_MAX 2048
#endif

/// @brief Largest real-SSA size emitted as chunked straight-line code (above it: interpreter).
#ifndef NT_CHUNK_MAX
#define NT_CHUNK_MAX 16384
#endif

/// @brief Straight-line slots per `noinline` block in the chunked tier (measured sweet spot:
///        matches unrolled speed within a block while bounding the optimiser's basic block).
#ifndef NT_CHUNK_BLOCK
#define NT_CHUNK_BLOCK 256
#endif

namespace numtracer::expr
{
  inline constexpr int kUnrollMax = NT_UNROLL_MAX;
  inline constexpr int kChunkMax = NT_CHUNK_MAX;
  inline constexpr int kChunkBlock = NT_CHUNK_BLOCK;
} // namespace numtracer::expr

// ---- portable compiler pragmas -------------------------------------------------------------
// Region guard that disables loop auto-vectorization for the enclosed function definitions.
// Sole use: working around a GCC 16 -O3 -march=native (AVX-512VL) miscompile of the complex-`Cx`
// std::vector fold loops in dense/dtensor.hpp (see there). Clang -O3 -march=native is NOT affected,
// so its branch is a documented no-op; any other compiler gets the default no-op. Bracket the
// affected definitions with `NT_BEGIN_NO_LOOP_VECTORIZE` ... `NT_END_NO_LOOP_VECTORIZE`.
#if defined(__GNUC__) && !defined(__clang__)
#define NT_BEGIN_NO_LOOP_VECTORIZE _Pragma("GCC push_options") _Pragma("GCC optimize(\"no-tree-loop-vectorize\")")
#define NT_END_NO_LOOP_VECTORIZE _Pragma("GCC pop_options")
#elif defined(__clang__)
#define NT_BEGIN_NO_LOOP_VECTORIZE // clang does not exhibit the GCC AVX-512VL fold bug
#define NT_END_NO_LOOP_VECTORIZE
#else
#define NT_BEGIN_NO_LOOP_VECTORIZE // unknown compiler: no-op
#define NT_END_NO_LOOP_VECTORIZE
#endif