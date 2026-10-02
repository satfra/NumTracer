/// @file core/config.hpp
/// @brief Portability shim: the exception-optional failure guard `NT_THROW`.
///
/// Nothing here changes *results* — only how a tripped guard is reported.
#pragma once

#include <cstddef>

// ---- exception-optional failure guard ------------------------------------------------------
// The library's internal misuse guards throw by default. The build-time net-builder generator
// TUs, however, are compiled with `-fno-exceptions` (see mathematica/CodegenBuild.m): emitting the
// exception-cleanup landing pads for their tens of thousands of destructible temporaries is what
// dominates their compile (docs/NUMTRACER_DESIGN_NOTES.md). Under `-fno-exceptions` a bare `throw` is ill-formed, so guards route through NT_THROW,
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
