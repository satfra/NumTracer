/// @file interpret.hpp
/// @brief Run a lowered program (@ref numtracer::GenProg) in-process — the arithmetic its emitted C++
///        performs, without compiling it. This is how a lowering is validated: interpret the program
///        and compare with evaluating the polynomial it came from.
///
/// Compare with a tolerance, not bitwise: the emitter fuses multiply-adds into `fma`, which rounds
/// differently from the separate `*` and `+` executed here.
#pragma once

#include "numtracer/codegen/gen.hpp" // GenProg, kRealProgram
#include "numtracer/codegen/real_cse.hpp"
#include "numtracer/core/cx.hpp"

#include <cstddef>
#include <vector>

namespace numtracer::inline network
{

  /// @brief The value of slot @p root of the instruction stream @p ins, reading variables from @p f
  ///        (`f[i]` = env slot i). A negative slot is a structural zero.
  inline double interpret(const std::vector<RInstr> &ins, int root, const double *f)
  {
    if (root < 0) return 0.0;
    std::vector<double> v(ins.size(), 0.0);
    auto val = [&](int r) { return r < 0 ? 0.0 : v[static_cast<std::size_t>(r)]; };
    for (std::size_t i = 0; i <= static_cast<std::size_t>(root); ++i) {
      const RInstr &in = ins[i];
      switch (in.op) {
      case RCONST: v[i] = in.value; break;
      case RVAR: v[i] = f[in.a]; break;
      case RADD: v[i] = val(in.a) + val(in.b); break;
      case RMUL: v[i] = val(in.a) * val(in.b); break;
      default: v[i] = -val(in.a); break;
      }
    }
    return val(root);
  }

  /// @brief The value of a lowered program at the env values @p f (from @ref Frame::fill_values).
  inline Cx interpret(const GenProg &prog, const double *f)
  {
    const double re = interpret(prog.ins, prog.root, f);
    const double im = prog.rootIm == kRealProgram ? 0.0 : interpret(prog.ins, prog.rootIm, f);
    return Cx{re, im};
  }

} // namespace numtracer::network
