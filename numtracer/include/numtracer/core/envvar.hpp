/// @file core/envvar.hpp
/// @brief The ONE truth test for every `NT_*` environment variable the build-time engine reads.
///
/// Everything here is generator-side configuration: the emitted kernel never reads an environment
/// variable, so nothing in this header is on a runtime hot path. Each reader caches its own value
/// in a function-local `static` at the call site — emission and contraction must be uniform across
/// a run, so a variable changed mid-process deliberately has no effect.
///
/// The rule: a flag is ON when the variable is set, non-empty, and not the single character "0".
/// Every reader must use it — a presence-only test reads `FOO=` and `FOO=0` as ON, a leading-char
/// test ignores `FOO=true` — and it must match `ntEnvFlag` in mathematica/DSL.m.
#pragma once

#include <cstdlib>
#include <cstring>

namespace numtracer
{

  /// @brief Read @p name as a boolean: set, non-empty, and not `"0"`.
  inline bool env_flag(const char *name)
  {
    const char *e = std::getenv(name);
    return e != nullptr && *e != '\0' && std::strcmp(e, "0") != 0;
  }

  /// @brief Read @p name as an integer, falling back to @p dflt when unset, empty or unparsable.
  ///
  /// Deliberately NOT `env_flag`-shaped: `0` is a legitimate value for several of these knobs
  /// (`NT_GEN_NOINLINE_MIN=0` puts every device trace function out of line), so emptiness — not the value — is
  /// what means "unset". Callers that need a range check apply it to the returned value; passing
  /// an out-of-range value is a caller error, not a parse error, and each call site says so.
  inline long env_int(const char *name, long dflt)
  {
    const char *e = std::getenv(name);
    if (e == nullptr || *e == '\0') return dflt;
    char *end = nullptr;
    const long v = std::strtol(e, &end, 10);
    return (end == e) ? dflt : v; // no digits consumed: treat as unset rather than as 0
  }

} // namespace numtracer
