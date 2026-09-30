#pragma once
/// @file precision.hpp
/// @brief Floating-point precision of EMITTED code.
///
/// The generator always computes in double; this only decides what the printed kernel computes in.
/// In single precision every emitted type is `float` and every literal carries an `f` suffix — a
/// bare double literal would silently promote the surrounding float arithmetic back to double.
/// Double precision leaves the emitted bytes exactly as they were before this switch existed.

#include "numtracer/core/config.hpp" // NT_THROW

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace numtracer::codegen
{
  enum class EmitPrecision { Double, Single };

  /// Process-wide, set once by the generator's `-p` flag before anything is emitted.
  inline EmitPrecision &emit_precision()
  {
    static EmitPrecision precision = EmitPrecision::Double;
    return precision;
  }

  inline bool emit_single() { return emit_precision() == EmitPrecision::Single; }

  /// The emitted real type.
  inline const char *emit_real_type() { return emit_single() ? "float" : "double"; }

  /// A float literal of @p v: rounded to float, printed with enough digits to round-trip it.
  inline std::string float_literal(double v)
  {
    const float fv = static_cast<float>(v);
    if (!std::isfinite(fv))
      NT_THROW(std::runtime_error, "float_literal: constant overflows float");
    std::ostringstream os;
    os.precision(std::numeric_limits<float>::max_digits10);
    os << fv;
    std::string s = os.str();
    if (s.find_first_of(".e") == std::string::npos) s += '.'; // "24f" is not a literal, "24.f" is
    return s + 'f';
  }
} // namespace numtracer::codegen
