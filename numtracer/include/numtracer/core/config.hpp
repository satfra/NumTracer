/// @file core/config.hpp
/// @brief Centralised compile-time build tunables and portability shims.
///
/// Gathered here so the trade-offs live in one place rather than being scattered through the
/// headers that use them. Nothing here changes *results* — only the build's cost profile or how
/// a tripped guard is reported.
#pragma once

#include <cstddef>

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