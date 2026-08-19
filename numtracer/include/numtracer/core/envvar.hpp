/// @file core/envvar.hpp
/// @brief The ONE truth test for every `NT_*` environment variable the build-time engine reads.
///
/// Everything here is generator-side configuration: the emitted kernel never reads an environment
/// variable, so nothing in this header is on a runtime hot path. Each reader caches its own value
/// in a function-local `static` at the call site — emission and contraction must be uniform across
/// a run, so a variable changed mid-process deliberately has no effect.
///
/// WHY THIS EXISTS. Four mutually incompatible spellings of "is this flag set" used to coexist in
/// the headers: `getenv(x) != nullptr` (presence-only), `e && e[0] == '1'` (exact leading char),
/// an open-coded copy of the rule below, and the rule itself, which was private to
/// `codegen/gen.hpp` and therefore unreachable from `numeric/`. Under the presence-only spelling
/// `FOO=` and `FOO=0` both read as TRUE, so *unsetting a flag by emptying it silently left it on*
/// — and under the leading-char spelling `FOO=true` silently did nothing. The Wolfram front end
/// solved the same problem with `ntEnvFlag` (mathematica/DSL.m) after `NT_NO_LABEL_CHECK=0`
/// disabled a correctness guard; this is the C++ half of that fix, and the two agree on the rule.
///
/// The rule: a flag is ON when the variable is set, non-empty, and not the single character "0".
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
  /// (`NT_GEN_SNAP_DIGITS=0` disables coefficient snapping), so emptiness — not the value — is
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
