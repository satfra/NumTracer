// nt_regulators.hpp — the regulators the generated kernels no longer carry.
//
// NumTracer emits a PLAIN kernel class: a flow whose dressing rules mention RB/RF/... emits
// unqualified calls to them and the consumer supplies the definitions. This header is the test
// suite's supply, and the single source of truth for the formulas — DiFfRG::ShimRegulator (shim.hpp)
// forwards here, so the NumTracer kernel and the still-templated FORM oracles it is graded
// against can never regulate differently.
//
// Definitions live at GLOBAL scope on purpose: unqualified lookup from a kernel member runs
// class -> kernel namespace -> global, so the same header serves kernels emitted into `DiFfRG`
// (the shim flows) and into `numtracer_kernels` (the standalone ones), with no per-namespace copy.
//
// The specific formulas are arbitrary but FIXED: the tests check that two evaluation PATHS agree as
// algebraic functions, so any consistent definitions used on both sides suffice. Exponential-type
// bosonic regulator, with nonzero analytic fermionic companions so the regularized quark propagator
// and its regulator-dot are genuinely exercised.
#pragma once

#include <cmath>

inline double RB(double k2, double p2) { return k2 * std::exp(-p2 / k2); }
inline double RBdot(double k2, double p2) { return 2.0 * (k2 + p2) * std::exp(-p2 / k2); }
inline double RF(double k2, double p2) { return std::sqrt(k2) * std::exp(-p2 / k2); }
inline double RFdot(double k2, double p2) { return (k2 + 2.0 * p2) / std::sqrt(k2) * std::exp(-p2 / k2); }
inline double dq2RB(double, double) { return 0.0; }
inline double dq2RF(double, double) { return 0.0; }
