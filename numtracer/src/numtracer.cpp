/// @file numtracer.cpp
/// @brief The one translation unit that *defines* the split engine/generator bodies (compiled mode).
///
/// Each engine header keeps its declarations + data-struct layouts + trivial builders always visible,
/// and wraps its heavy function definitions in `#if NUMTRACER_DEFINE_BODIES`. Defining
/// `NUMTRACER_COMPILING_LIB` here turns those regions on for this TU only, so the bodies are compiled
/// and optimised a single time and linked into `NumTracer::NumTracer`. In a normal compiled consumer
/// TU the regions are skipped (declarations only → link to this object). In the header-only build
/// (`NUMTRACER_HEADER_ONLY=1`) the bodies are `inline` in the headers and this TU is effectively empty.
/// See core/export.hpp.
#define NUMTRACER_COMPILING_LIB 1

// Pulls the whole engine: numeric contraction + driver, SU(N) folding, lowering/CSE, polynomial
// arithmetic. With NUMTRACER_COMPILING_LIB set, every split header emits its definitions here once.
#include "numtracer/numeric/numeric_driver.hpp"
#include "numtracer/network/sun_net.hpp"
#include "numtracer/codegen/gen.hpp"

