#!/usr/bin/env python3
"""host_test.py -- build a host-only test of a flow rewritten by emit_tables.py.

For random nodes, the table entries from tab::coeffs() are fed to the device reconstruction r<trN>() and compared
with the emitted float32 trace and a long-double reference (the same emitted code, evaluated in tab::), over
S0 in [0.005, 250] and k in [0.034, 1000].  Metric per (trace, S0, k): max over nodes |T - T_ref| / max |T_ref|;
reported as the max over all cells and over the bulk 1e-2 <= S0/k <= 1e2.

Usage: host_test.py <rewritten_flow_dir> <out.cc>   then   g++-14 -O2 -std=c++20 out.cc && ./a.out
"""
import re
import sys
from pathlib import Path

flow = Path(sys.argv[1]).resolve()
khh = (flow / "kernels.hh").read_text()
ns = re.search(r"namespace DiFfRG \{ namespace (\w+) \{", khh).group(1)
tabbed = re.findall(r"KOKKOS_INLINE_FUNCTION float r(tr\d+)\(", khh)
args = [a.split()[-1] for a in re.search(r"void fill\(float \*f, ([^)]*)\)", khh).group(1).split(",")]
shape = [a.split()[-1] for a in re.search(r"inline void coeffs\(float \*out, double qhat, ([^)]*)\)", khh).group(1).split(",")]
scale = "S0" if "S0" in args else "p"
draw = {"cos1": "-0.99 + 1.98 * U(g)", "cos2": "-0.99 + 1.98 * U(g)", "S1": "0.99 * U(g)", "SPhi": "6.28 * U(g)",
        "phi": "6.28 * U(g)"}
calls = []
for i, t in enumerate(tabbed):
    calls.append(f"      agg({i}, cell, tab::{t}(fr).re, {t}(fenv), r{t}(row.data(), 1, fenv, s, r, kle));")
src = f"""#define KOKKOS_INLINE_FUNCTION inline
#define KOKKOS_FUNCTION
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>
#include "{flow}/kernels.hh"
using namespace DiFfRG::{ns};
constexpr int NT = {len(tabbed)};
static const double S0s[] = {{0.005, 0.05, 0.5, 5., 50., 250.}}, ks[] = {{0.034, 0.3, 3., 30., 300., 1000.}};
static long double mref[NT][36], m32[NT][36], mtab[NT][36];
static void agg(int i, int c, long double ref, double a, double b) {{
  mref[i][c] = std::fmax(mref[i][c], std::fabs(ref));
  m32[i][c] = std::fmax(m32[i][c], std::fabs(a - ref));
  mtab[i][c] = std::fmax(mtab[i][c], std::fabs(b - ref));
}}
int main() {{
  std::mt19937 g(5); std::uniform_real_distribution<double> U(0, 1);
  std::vector<float> row(8192);
  for (int node = 0; node < 200; ++node) {{
    const double qhat = 0.01 + 2.5 * U(g);
    {"".join(f"const double {a} = {draw[a]}; " for a in shape)}
    tab::coeffs(row.data(), qhat, {", ".join(shape)});
    for (int is = 0; is < 6; ++is) for (int ik = 0; ik < 6; ++ik) {{
      const double {scale} = S0s[is], k = ks[ik], l1 = k * qhat; const int cell = is * 6 + ik;
      float fenv[nenv]; fill(fenv, {", ".join(f"float({a})" for a in args)});
      tab::Cx fr[nenv]; tab::fill(fr, {", ".join("tab::Cx(l1)" if a == "l1" else f"tab::R({a})" for a in args)});
      const float s = k > {scale} ? k : {scale}; const float r = (k > {scale} ? {scale} : k) / s; const bool kle = k <= {scale};
{chr(10).join(calls)}
    }}
  }}
  double all32 = 0, alltab = 0, bulk32 = 0, bulktab = 0;
  for (int i = 0; i < NT; ++i) for (int c = 0; c < 36; ++c) if (mref[i][c] > 0) {{
    const double ratio = S0s[c / 6] / ks[c % 6];
    const double e32 = m32[i][c] / mref[i][c], et = mtab[i][c] / mref[i][c];
    all32 = std::fmax(all32, e32); alltab = std::fmax(alltab, et);
    if (ratio >= 1e-2 && ratio <= 1e2) {{ bulk32 = std::fmax(bulk32, e32); bulktab = std::fmax(bulktab, et); }}
  }}
  std::printf("{flow.name:14s} all cells: inline32 %.1e table %.1e | bulk 1e-2..1e2: inline32 %.1e table %.1e\\n",
              all32, alltab, bulk32, bulktab);
}}
"""
Path(sys.argv[2]).write_text(src)
