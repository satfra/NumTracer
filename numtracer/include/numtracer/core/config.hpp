/// @file core/config.hpp
/// @brief Centralised compile-time build tunables and portable compiler-pragma shims.
///
/// Gathered here so the trade-offs live in one place rather than being scattered through the
/// headers that use them. Nothing here changes *results* — only the build's cost profile or a
/// compiler-specific workaround.
#pragma once

#include <cstddef>

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

// ---- exception-optional failure guard ------------------------------------------------------
// The library's internal misuse guards throw by default. The build-time net-builder generator
// TUs, however, are compiled with `-fno-exceptions` (see mathematica/CodegenBuild.m): emitting the
// exception-cleanup landing pads for their tens of thousands of destructible temporaries is what
// dominates their -O0 compile — turning exceptions off cut a representative unit from 15.3 s to
// 1.3 s. Under `-fno-exceptions` a bare `throw` is ill-formed, so guards route through NT_THROW,
// which keeps the exact same exception (type + message) when exceptions are enabled and degrades
// to a loud abort() when they are not. A tripped guard is always a bug, never a recoverable
// condition — correct runs never reach it — so results are identical either way. `exc` is the
// std::exception subclass (the caller includes <stdexcept>); `msg` is a `const char*`.
#if defined(__cpp_exceptions) || defined(__EXCEPTIONS)
#define NT_THROW(exc, msg) throw exc(msg)
#else
#include <cstdio>
#include <cstdlib>
#define NT_THROW(exc, msg) (std::fprintf(stderr, "numtracer fatal: %s\n", msg), std::abort())
#endif