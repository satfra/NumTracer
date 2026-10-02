/// @file numtracer.hpp
/// @brief Convenience umbrella header: pulls in the whole public NumTracer API so a consumer
///        can `#include <numtracer.hpp>` instead of the individual subsystem headers.
///
/// This bundles the numeric contraction engine (the hand-coding path: declare a @ref numtracer::Frame,
/// build a `DiracChain` + `LorentzNet`, contract with `Frame::trace`, evaluate the resulting `Poly`),
/// the SU(N) fold (@ref numtracer::SUN), and the build-time codegen/lowering helpers. Every header below is `#pragma once`-guarded, so
/// including this alongside any individual header is harmless.
///
/// It does NOT change how the engine ships: the compiled-vs-header-only choice is still governed by
/// `numtracer/core/export.hpp` (`NUMTRACER_HEADER_ONLY`) — a binary must pick one mode across its
/// TUs regardless of whether it includes this umbrella or the individual headers.
#pragma once

// Every include below is deliberately re-exported (this header's whole purpose); the trailing
// `IWYU pragma: export` markers say so, so include-what-you-use / clangd do not flag them "unused".

// Core: the constexpr complex NTTP, the compile-time constant carrier, the dense complex matrix,
// and the compiled-vs-header-only export gate.
#include "numtracer/core/cx.hpp"      // IWYU pragma: export
#include "numtracer/core/lit.hpp"     // IWYU pragma: export
#include "numtracer/core/cmat.hpp"    // IWYU pragma: export
#include "numtracer/core/export.hpp"  // IWYU pragma: export

// Typed-out tables (cross-checked against the runtime oracles): Euclidean Weyl gamma matrices and
// the SU(2)/SU(3) structure/generator data.
#include "numtracer/dirac/dirac_data.hpp" // IWYU pragma: export
#include "numtracer/sun/sun_data.hpp"      // IWYU pragma: export

// Tensor network: the generic network, the Dirac-chain tokens, and the numeric SU(N) fold.
#include "numtracer/network/network.hpp" // IWYU pragma: export
#include "numtracer/network/dirac.hpp"   // IWYU pragma: export
#include "numtracer/network/sun_net.hpp" // IWYU pragma: export

// Numeric contraction engine: the polynomial types, the 4x4 chiral spinor trace, the contractor,
// and the all-diagram driver.
#include "numtracer/numeric/mpoly.hpp"            // IWYU pragma: export
#include "numtracer/numeric/dpoly.hpp"            // IWYU pragma: export
#include "numtracer/numeric/spinor_mat.hpp"       // IWYU pragma: export
#include "numtracer/numeric/numeric_contract.hpp" // IWYU pragma: export
#include "numtracer/numeric/numeric_driver.hpp"   // IWYU pragma: export
#include "numtracer/numeric/frame.hpp"            // IWYU pragma: export

// Codegen: the real-SSA CSE, the Poly lowering, the emitter, and the emitted-kernel support layer
// (numtracer::complex / numtracer::compute) that a generated kernel compiles against.
#include "numtracer/codegen/real_cse.hpp" // IWYU pragma: export
#include "numtracer/codegen/lower.hpp"    // IWYU pragma: export
#include "numtracer/codegen/gen.hpp"      // IWYU pragma: export
#include "numtracer/codegen/runtime.hpp"  // IWYU pragma: export
#include "numtracer/codegen/interpret.hpp" // IWYU pragma: export
