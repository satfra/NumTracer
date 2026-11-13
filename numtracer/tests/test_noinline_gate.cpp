// NumTracer — the device `__noinline__` size gate, tested through its only observable: the
// decorator `emit_cpp` actually prints.
//
// WHY THIS EXISTS. The gate (gen.hpp's `edetail::eff_decor`) was dead in production for months, in
// two successive ways, and neither failure produced a diagnostic — it emitted a perfectly valid
// kernel that merely happened to be all-inline:
//
//   1. it sniffed the decorator for a literal `__device__`, which `ntKokkosDecor` has already
//      rewritten to the Kokkos macros by then, so the test never matched;
//   2. once that was replaced by an explicit `NT_GEN_DEVICE`, the variable was set only as a shell
//      prefix on the ONLINE Wolfram `Run[]` path — while every production flow is generated
//      OFFLINE, from a `cmake -P` step that inherits none of that environment.
//
// A gate whose failure mode is "silently does nothing" needs a test that asserts it DOES something,
// which is what the `fires` case below is for. The rest pin the three ways to ask for no gating,
// one of which is a genuine trap: the threshold test is `nInstr > N`, so `NT_GEN_NOINLINE_MIN=0`
// out-of-lines EVERYTHING. The escape hatch back to all-inline emission is `off`, not `0`.
//
// The environment is read once per process and cached (deliberately — emission must be consistent
// across every function in a run), so each mode needs its own process. Hence: one binary, the mode
// named on argv, and ctest registering one case per mode with its own ENVIRONMENT.
#include "numtracer/codegen/gen.hpp"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>

using namespace numtracer;
using namespace numtracer::network;

namespace
{

  /// A straight-line program of @p nInstr instructions, rooted at the last one. The arithmetic is
  /// irrelevant — the gate keys on `p.ins.size()` alone — but it must be a real chain so the
  /// emitter's liveness plan keeps every slot and prints one statement each.
  GenProg chain(std::size_t nInstr)
  {
    GenProg p;
    p.ins.push_back({RVAR, 0, -1, 0.0});
    for (std::size_t i = 1; i < nInstr; ++i)
      p.ins.push_back({RMUL, static_cast<int>(i - 1), 0, 0.0});
    p.root = static_cast<int>(p.ins.size()) - 1;
    return p;
  }

  /// The decorator `emit_cpp` printed for a program of @p nInstr instructions — i.e. everything up
  /// to the return type on the emitted signature line.
  std::string decorated(std::size_t nInstr)
  {
    std::ostringstream os;
    emit_cpp(os, chain(nInstr), "tr0", "static KOKKOS_INLINE_FUNCTION");
    const std::string s = os.str();
    const std::size_t at = s.find(" double tr0(");
    return at == std::string::npos ? std::string("<no signature emitted>") : s.substr(0, at);
  }

  int fail = 0;

  void ok(const char *what, bool cond)
  {
    std::printf("%-58s %s\n", what, cond ? "ok" : "FAIL");
    if (!cond) ++fail;
  }

  bool outOfLined(const std::string &decor)
  {
    return decor.find("__attribute__((noinline))") != std::string::npos;
  }

} // namespace

int main(int argc, char **argv)
{
  const std::string mode = argc > 1 ? argv[1] : "host";

  // Sizes straddling the default threshold of 500. 600 is over it, 100 well under.
  const std::string big = decorated(600), small = decorated(100);
  std::printf("mode=%s\n  600-instr: %s\n  100-instr: %s\n\n", mode.c_str(), big.c_str(), small.c_str());

  if (mode == "host") {
    // No NT_GEN_DEVICE: the host has no register cliff and its emission must stay byte-identical to
    // what it always was, whatever the size.
    ok("host: large function stays inline", !outOfLined(big));
    ok("host: small function stays inline", !outOfLined(small));
    ok("host: keeps the Kokkos INLINE spelling", big.find("KOKKOS_INLINE_FUNCTION") != std::string::npos);
    ok("host: a trace is a template (compiled only where called)", big.rfind(kLazyPrefix, 0) == 0);

  } else if (mode == "fires") {
    // NT_GEN_DEVICE=1, default threshold. THE regression this file exists for: if the gate ever
    // goes dead again, this is what catches it.
    ok("device: large function is out-of-lined", outOfLined(big));
    ok("device: small function stays inline", !outOfLined(small));
    // KOKKOS_INLINE_FUNCTION expands to __forceinline__, which would contradict the attribute; the
    // gate must swap in the plain spelling rather than stack the two.
    ok("device: swaps in KOKKOS_FUNCTION", big.find("KOKKOS_FUNCTION") != std::string::npos);
    ok("device: drops the INLINE spelling", big.find("KOKKOS_INLINE_FUNCTION") == std::string::npos);

  } else if (mode == "min0") {
    // The trap. `nInstr > N` means 0 gates EVERYTHING with any instruction at all — it is not the
    // way to turn the gate off, and reading it as such would silently out-of-line every flow.
    ok("min=0: large function is out-of-lined", outOfLined(big));
    ok("min=0: SMALL function is ALSO out-of-lined", outOfLined(small));

  } else if (mode == "off") {
    // The actual escape hatch.
    ok("off: large function stays inline", !outOfLined(big));
    ok("off: small function stays inline", !outOfLined(small));
    ok("off: keeps the Kokkos INLINE spelling", big.find("KOKKOS_INLINE_FUNCTION") != std::string::npos);

  } else if (mode == "forced") {
    // NT_GEN_NOINLINE_TRACES forces out-of-line everywhere, host included, ignoring the threshold.
    ok("forced: large function is out-of-lined", outOfLined(big));
    ok("forced: small function is out-of-lined too", outOfLined(small));

  } else {
    std::printf("unknown mode '%s'\n", mode.c_str());
    return 2;
  }

  std::printf("\n%s (%d failure%s)\n", fail ? "TESTS FAILED" : "ALL TESTS PASSED", fail,
              fail == 1 ? "" : "s");
  return fail ? 1 : 0;
}
