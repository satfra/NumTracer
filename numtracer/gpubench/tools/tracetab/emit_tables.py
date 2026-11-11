#!/usr/bin/env python3
"""emit_tables.py -- PROTOTYPE: rewrite an emitted flow so its traces come from a per-run coefficient table.

Every trace T is homogeneous in (l1, S) (S = the dimensionful external coordinate).  With l1 = k*qhat and the inv
slots cleared (Q_j = 1/inv_j, mass dimension w_j, max power m_j in the trace),
    N = T * prod Q_j^m_j = sum_{n=n0..n1} c_n k^n S^(D-n),   D = d + sum_j w_j m_j,
where c_n depends only on (qhat, angles, shape) -- not on k or S.  At run time, with s = max(k, S), r = min/max:
    T = s^d * prod_j (s^w_j inv_j)^m_j * r^(n0 or D-n1) * Horner_r(c)   (coefficient order reversed for k > S).
Only the live range n0..n1 is stored: coefficients that vanish identically would otherwise come back as rounding
noise, amplified by (k/S)^n at the extremes.

The coefficients are computed on the host when DiFfRG builds the table (diffrg_tab.patch): N is sampled on a circle
in the complex l1 plane (away from the zeros of the Q_j), in long double, and an inverse DFT gives its coefficients.
The structure (m_j, w_j, d, n0, n1) comes from probing the emitted code with analyse.P.  Since the emitted float32
code carries float-rounded constants, heavily cancelling traces are only as exact as those constants; a production
version must take the coefficients from NumTracer's exact polynomials instead.

Usage: emit_tables.py <flow_dir> <out_dir> [--balance B] [--scale NAME]
  A trace is tabulated when its inline flop count exceeds B * 4 * (n1-n0+1) (bytes read per node); default B = 1.
"""
from __future__ import annotations

import argparse
import math
import random
import re
from collections import defaultdict
from pathlib import Path

import analyse
from analyse import P
from validate import compile_fill, compile_trace, fma64

NUM = analyse.num


def complexify_expr(e: str) -> str:
    e = NUM.sub(r"\1", e)
    return e.replace("powr<", "powrX<").replace("fma(", "fmaX(")


def structure(slots, traces, args, fill, radial, scale, rnd):
    """Per trace: (name, flop, d, {inv slot: (m_j, w_j)}, D, n0, n1)."""
    rad, inv, other = analyse.classify(slots)
    P.nv = 1 + len(inv)
    sidx = {j: 1 + n for n, j in enumerate(inv)}
    shape_names = [a for a in args if a not in (radial, scale)]

    def point():
        d = {}
        for a in shape_names:
            d[a] = rnd.uniform(0.0, 0.99) if a == "S1" else rnd.uniform(-0.95, 0.95) if a.startswith("cos") \
                else rnd.uniform(0.1, 6.2)
        return d

    # mass dimension of each inv slot
    base = point()
    fa = fill(**{**base, radial: 0.7, scale: 1.0})
    fb = fill(**{**base, radial: 1.4, scale: 2.0})
    w = {j: round(math.log2(fa[j] / fb[j])) for j in inv}
    assert all(v in (0, 2) for v in w.values()), w

    out = []
    for name, body, ret in traces:
        tr = compile_trace(body, ret)
        m = defaultdict(int)
        polys = []
        for _ in range(3):
            pt = point()
            vals = fill(**{**pt, radial: 1.0, scale: 1.0})
            env = {j: (P.var(0, 1.0) if j == rad[0] else P.var(sidx[j]) if j in sidx else vals[j])
                   for j in range(len(vals))}
            poly = tr(env, fma64)
            poly = poly if isinstance(poly, P) else P.const(float(poly))
            # Q_j(l1) = 1/inv_j at S = 1, a quadratic in l1 (constant for the l1-independent slots)
            q = [fill(**{**pt, radial: x, scale: 1.0}) for x in (0.0, 1.0, 2.0)]
            Q = {}
            for j in inv:
                y = [1.0 / q[i][j] for i in range(3)]
                c2 = (y[2] - 2 * y[1] + y[0]) / 2
                Q[j] = [y[0], y[1] - y[0] - c2, c2]
            polys.append((poly, Q))
            for mono in poly.live():
                for j in inv:
                    m[j] = max(m[j], mono[sidx[j]])
        # the k-power range of the cleared numerator N = T prod Q_j^m_j: exact zeros at either end must be dropped,
        # since they are amplified by (k/S)^n at the extremes
        n0, n1 = 10 ** 9, -1
        for poly, Q in polys:
            Nc = defaultdict(lambda: [0.0, 0.0])
            for mono, (c, am) in poly.t.items():
                terms = {mono[0]: [c, am]}
                for j in inv:
                    for _ in range(m[j] - mono[sidx[j]]):
                        nt = defaultdict(lambda: [0.0, 0.0])
                        for dd, (v, va) in terms.items():
                            for e, qe in enumerate(Q[j]):
                                nt[dd + e][0] += v * qe
                                nt[dd + e][1] += va * abs(qe)
                        terms = nt
                for dd, (v, va) in terms.items():
                    Nc[dd][0] += v
                    Nc[dd][1] += va
            live_n = [dd for dd, (v, va) in Nc.items() if va > 0 and abs(v) > 1e-6 * va]
            if live_n:
                n0, n1 = min(n0, min(live_n)), max(n1, max(live_n))
        votes = defaultdict(int)
        for _ in range(12):
            pt = point()
            t1 = tr(fill(**{**pt, radial: 0.8, scale: 1.0}), fma64)
            t2 = tr(fill(**{**pt, radial: 1.6, scale: 2.0}), fma64)
            if t1 != 0 and t2 != 0:
                votes[round(math.log2(abs(t2 / t1)))] += 1
        d = max(votes, key=votes.get) if votes else 0
        flop = count_flop(body)
        mj = {j: (m[j], w[j]) for j in inv if m[j] > 0}
        D = d + sum(mm * ww for mm, ww in mj.values())
        if n1 < 0:  # identically zero
            n0 = n1 = 0
        assert 0 <= n0 <= n1 <= D, (name, n0, n1, D)
        out.append((name, flop, d, mj, D, n0, n1))
    return out


def count_flop(body):
    n = 0
    for line in body.split("\n"):
        line = line.strip().replace("[[maybe_unused]] ", "")
        if not line.startswith("const "):
            continue
        rhs = line.split(" = ", 1)[1]
        if rhs.startswith("fma("):
            n += 2
        elif re.search(r"[*+]|(?<!\()-", rhs) and not re.fullmatch(r"\(?-?[\d.ef+-]+\)?;", rhs):
            n += 1
    return n


def find_function(text, header_re):
    """(start, end) of a member function whose signature line matches header_re, including its template line."""
    m = re.search(header_re, text)
    if not m:
        return None
    start = text.rfind("\n", 0, m.start()) + 1
    prev = text.rfind("\n", 0, start - 1) + 1
    if text[prev:start].strip().startswith("template"):
        start = prev
    i = text.index("{", m.end() - 1)
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return start, j + 1
    raise ValueError("unbalanced")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("flow_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--balance", type=float, default=1.0)
    ap.add_argument("--scale", default=None)
    a = ap.parse_args()

    khh = (a.flow_dir / "kernels.hh").read_text()
    kh = (a.flow_dir / "kernel.hh").read_text()
    ns = re.search(r"namespace DiFfRG \{ namespace (\w+) \{", khh).group(1)
    slots, traces = analyse.parse(str(a.flow_dir / "kernels.hh"))
    rad, inv, other = analyse.classify(slots)
    args, fill = compile_fill(khh)
    radial = slots[rad[0]]
    scale = a.scale or ("S0" if "S0" in args else "p")
    assert scale in args, (scale, args)
    rnd = random.Random(17)
    st = structure(slots, traces, args, fill, radial, scale, rnd)

    tab, off = [], 0
    for (name, flop, d, mj, D, n0, n1) in st:
        if flop > a.balance * 4 * (n1 - n0 + 1):
            tab.append((name, flop, d, mj, D, n0, n1, off))
            off += n1 - n0 + 1
    E = off
    n1max = max(t[6] for t in tab)
    print(f"{a.flow_dir.name}: {len(tab)}/{len(st)} traces tabulated, {E} entries per node, max k-power {n1max}, "
          f"inline flop {sum(t[1] for t in st)} -> kept inline {sum(t[1] for t in st) - sum(t[1] for t in tab)}")

    # ---- kernels.hh: host coefficient code + device reconstruction ----
    fillm = re.search(r"void fill\(\w+ \*f, ([^)]*)\) \{(.*?)\n\}", khh, re.S)
    fbody = []
    for i, rhs in re.findall(r"f\[(\d+)\] = (.*?);\s*$", fillm.group(2), re.M):
        fbody.append(f"  f[{i}] = Cx({complexify_expr(rhs)});")
    cargs = ", ".join(("Cx " if n == radial else "R ") + n for n in args)
    tbodies = {n: (b, r) for n, b, r in traces}
    host = [f"\n// ---- trace tables (PROTOTYPE, emit_tables.py) ----",
            "#include <atomic>", "#include <cstdio>",
            f"namespace DiFfRG {{ namespace {ns} {{ namespace tab {{",
            "using R = long double; // extra digits: the emitted monomial form cancels heavily for some traces",
            "struct Cx { R re = 0, im = 0; Cx() = default; Cx(R r, R i = 0) : re(r), im(i) {} };",
            "inline Cx operator+(Cx a, Cx b) { return {a.re + b.re, a.im + b.im}; }",
            "inline Cx operator-(Cx a, Cx b) { return {a.re - b.re, a.im - b.im}; }",
            "inline Cx operator-(Cx a) { return {-a.re, -a.im}; }",
            "inline Cx operator*(Cx a, Cx b) { return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}; }",
            "inline Cx operator/(Cx a, Cx b) { const R n = b.re * b.re + b.im * b.im; "
            "return {(a.re * b.re + a.im * b.im) / n, (a.im * b.re - a.re * b.im) / n}; }",
            "inline Cx fmaX(Cx a, Cx b, Cx c) { return a * b + c; }",
            "inline R fmaX(R a, R b, R c) { return a * b + c; }",
            "template <int N, class T> inline T powrX(T x) { T r = T(1.); for (int i = 0; i < (N < 0 ? -N : N); ++i) r = r * x; "
            "if constexpr (N < 0) return T(1.) / r; else return r; }",
            "using std::sqrt; using std::cos; using std::sin;",
            f"inline void fill(Cx *f, {cargs}) {{", *fbody, "}"]
    for (name, *_rest) in tab:
        b, r = tbodies[name]
        host.append(f"inline Cx {name}(const Cx *f) {{")
        for line in b.split("\n"):
            line = line.strip().replace("[[maybe_unused]] ", "")
            if not line.startswith("const "):
                continue
            lhs, rhs = line.split(" ", 2)[2].rstrip(";").split(" = ", 1)
            host.append(f"  const Cx {lhs} = {complexify_expr(rhs)};")
        ret_expr = NUM.sub(r"\1", r.rstrip(";"))
        host.append(f"  return Cx({ret_expr});")
        host.append("}")
    # inv slots whose denominator depends on the loop momentum: their zeros constrain the sampling radius
    lin = [j for j in inv if any(j in t[3] for t in tab) and re.search(r"\b" + radial + r"\b", slots[j])]
    M = n1max + 3
    shape_args = [n for n in args if n not in (radial, scale)]
    host += [
        "inline std::atomic<double> worst_alias{0.}, worst_check{0.};",
        "inline void atomic_max(std::atomic<double> &a, double v) { double c = a.load(); while (v > c && !a.compare_exchange_weak(c, v)) {} }",
        "struct Report { ~Report() { std::fprintf(stderr, \"[tab] " + ns + ": worst aliasing %.2e, worst self-check %.2e\\n\", worst_alias.load(), worst_check.load()); } };",
        "inline Report report;",
        f"/// The {E} table entries of a node: k-coefficients of every tabulated trace at S = 1.",
        f"inline void coeffs(float *out, double qhat, {', '.join('double ' + n for n in shape_args)}) {{",
        "  constexpr R pi = 3.14159265358979323846264338327950288L;",
        f"  constexpr int M = {M};",
        f"  auto F = [&](Cx *f, Cx l) {{ fill(f, {', '.join('l' if n == radial else '1.0' if n == scale else n for n in args)}); }};",
        f"  Cx f[{len(slots)}];",
        "  // radius in the l1 plane: as far as possible from the zeros of the l1-dependent Q_j (|zero| = sqrt(c0/c2))",
        "  R Rz[" + str(max(1, len(lin))) + "]; int nR = 0;",
        "  {",
        "    Cx f0[" + str(len(slots)) + "], f1[" + str(len(slots)) + "], fm[" + str(len(slots)) + "];",
        "    F(f0, Cx(0.)); F(f1, Cx(1.)); F(fm, Cx(-1.));",
    ]
    for j in lin:
        host += [f"    {{ const R q0 = (Cx(1.) / f0[{j}]).re, q1 = (Cx(1.) / f1[{j}]).re, qm = (Cx(1.) / fm[{j}]).re;",
                 "      const R c2 = 0.5 * (q1 + qm) - q0; Rz[nR++] = std::sqrt(std::fabs(q0 / c2)); }"]
    host += [
        "  }",
        "  std::sort(Rz, Rz + nR);",
        "  R rho = 1., best = -1.;",
        "  auto score = [&](R r) { R s = 1e300; for (int i = 0; i < nR; ++i) s = std::fmin(s, std::fabs(std::log(r / Rz[i]))); return s; };",
        "  if (nR == 0) rho = 1.; else {",
        f"    R cand[{len(lin) + 2}]; int nc = 0; cand[nc++] = 0.5 * Rz[0]; cand[nc++] = 2. * Rz[nR - 1];",
        "    for (int i = 0; i + 1 < nR; ++i) cand[nc++] = std::sqrt(Rz[i] * Rz[i + 1]);",
        "    for (int i = 0; i < nc; ++i) if (score(cand[i]) > best) { best = score(cand[i]); rho = cand[i]; }",
        "  }",
        f"  Cx N[{len(tab)}][M];",
        "  for (int s = 0; s < M; ++s) {",
        "    const Cx z = Cx(rho * std::cos(2 * pi * s / M), rho * std::sin(2 * pi * s / M));",
        "    F(f, z);",
    ]
    for i, (name, flop, d, mj, D, n0, n1, o) in enumerate(tab):
        fac = " * ".join(f"powrX<{mm}>(Cx(1.) / f[{j}])" for j, (mm, ww) in mj.items()) or "Cx(1.)"
        host.append(f"    N[{i}][s] = {name}(f) * {fac};")
    host += ["  }"]
    # inverse DFT, aliasing check, self-check at two real points
    host += [
        "  Cx fr[" + str(len(slots)) + "];",
        "  for (int t = 0; t < " + str(len(tab)) + "; ++t) {",
        "    static constexpr int N0[] = {" + ", ".join(str(t[5]) for t in tab) + "};",
        "    static constexpr int N1[] = {" + ", ".join(str(t[6]) for t in tab) + "};",
        "    static constexpr int Os[] = {" + ", ".join(str(t[7]) for t in tab) + "};",
        "    R a[M], amax = 0.;",
        "    for (int n = 0; n < M; ++n) {",
        "      Cx acc{};",
        "      for (int s = 0; s < M; ++s) acc = acc + N[t][s] * Cx(std::cos(2 * pi * s * n / M), -std::sin(2 * pi * s * n / M));",
        "      a[n] = acc.re / M * std::pow(rho, -n);",
        "      amax = std::fmax(amax, std::fabs(acc.re) / M);",
        "    }",
        "    R alias = 0.;",
        "    for (int n = 0; n < M; ++n) if (n < N0[t] || n > N1[t]) alias = std::fmax(alias, std::fabs(a[n]) * std::pow(rho, n));",
        "    if (amax > 0) atomic_max(worst_alias, double(alias / amax));",
        "    for (int n = N0[t]; n <= N1[t]; ++n) out[Os[t] + n - N0[t]] = float(a[n] * std::pow(qhat, n));",
        "  }",
    ]
    # self check at l1 = qhat * kk, S = 1 (s = 1)
    host += ["  for (R kk : {0.3L, 1.0L}) {", "    F(fr, Cx(qhat * kk));"]
    for i, (name, flop, d, mj, D, n0, n1, o) in enumerate(tab):
        fac = " * ".join(f"std::pow(fr[{j}].re, {mm})" for j, (mm, ww) in mj.items()) or "1."
        host += [f"    {{ const R dir = {name}(fr).re; R h = 0., mag = 0.;",
                 f"      for (int n = {n1}; n >= {n0}; --n) {{ h = h * kk + out[{o} + n - {n0}]; mag = mag * kk + std::fabs(out[{o} + n - {n0}]); }}",
                 f"      h *= std::pow(kk, {n0}); mag *= std::pow(kk, {n0});",
                 f"      const R fac = {fac};",
                 "      if (mag * std::fabs(fac) > 0) atomic_max(worst_check, double(std::fabs(h * fac - dir) / (mag * std::fabs(fac)))); }"]
    host += ["  }", "}", "}}} // namespace tab"]

    dev = [f"namespace DiFfRG {{ namespace {ns} {{"]
    for (name, flop, d, mj, D, n0, n1, o) in tab:
        dev.append(f"static KOKKOS_INLINE_FUNCTION float r{name}(const float *row, const size_t st, const float *f, "
                   "const float s, const float r, const bool kle) {")
        dev.append(f"  // N = sum_(n={n0}..{n1}) c_n k^n S^({D}-n); row holds c_{n0}..c_{n1}")
        dev.append("  float h;")
        dev.append("  if (kle) {")
        dev.append(f"    h = row[{o + n1 - n0} * st];")
        for n in range(n1 - 1, n0 - 1, -1):
            dev.append(f"    h = fma(h, r, row[{o + n - n0} * st]);")
        dev.append(f"    h *= powr<{n0}>(r);" if n0 else "")
        dev.append("  } else {")
        dev.append(f"    h = row[{o} * st];")
        for n in range(n0 + 1, n1 + 1):
            dev.append(f"    h = fma(h, r, row[{o + n - n0} * st]);")
        dev.append(f"    h *= powr<{D - n1}>(r);" if D - n1 else "")
        dev.append("  }")
        fac = [f"powr<{d}>(s)"] if d else []
        for j, (mm, ww) in mj.items():
            base = f"(s * s * f[{j}])" if ww == 2 else f"f[{j}]"
            fac.append(f"powr<{mm}>({base})")
        dev.append(f"  return h{''.join(' * ' + x for x in fac)};")
        dev.append("}")
    dev.append("}}")

    new_khh = khh.replace("#include <cmath>", "#include <cmath>\n#include <algorithm>", 1)
    new_khh += "\n".join(host) + "\n" + "\n".join(dev) + "\n"

    # ---- kernel.hh ----
    tabbed = {t[0] for t in tab}
    sig = re.search(r"static KOKKOS_INLINE_FUNCTION auto kernel\((.*?)\)\s*\n", kh).group(1)
    params = [p.strip() for p in sig.split(",")]
    fl = [p for p in params if p.startswith("const float&")]
    names = [p.split("&")[1].strip() for p in fl]
    assert "k" in names, names
    names = names[: names.index("k") + 1]
    ext = [n for n in names[3:-1]]
    mask = ", ".join("false" if n == scale else "true" for n in ext)
    tab_fn = (f"    static constexpr bool has_tab = true;\n"
              f"    static constexpr size_t tab_size = {E};\n"
              f"    static constexpr std::array<bool, {len(ext)}> tab_ext_mask{{{mask}}};\n"
              f"    template <typename... Rest_>\n"
              f"    static void tab(float *out_, {', '.join('const float& ' + n for n in names)}, const Rest_ &...)\n"
              f"    {{ DiFfRG::{ns}::tab::coeffs(out_, double({radial}) / double(k), "
              f"{', '.join(n for n in shape_args)}); }}\n")

    def tabify(fn_text, kind):
        if kind == "kernel":
            fn_text = fn_text.replace("auto kernel(", "auto kernel_tab(const float *row_, const size_t stride_, ", 1)
        else:
            fn_text = fn_text.replace("auto post(const TP_ &P_, ", "auto post_tab(const TP_ &P_, const float *row_, const size_t stride_, ", 1)
        fn_text = re.sub(r"(\n\s*DiFfRG::" + ns + r"::fill\(fenv[^;]*;)",
                         r"\1\n      const float _ts = k > " + scale + " ? k : " + scale + "; const float _tr = (k > " + scale +
                         " ? " + scale + " : k) / _ts; const bool _tkle = k <= " + scale + ";", fn_text, count=1)
        for t in tabbed:
            fn_text = fn_text.replace(f"DiFfRG::{ns}::{t}(fenv)", f"DiFfRG::{ns}::r{t}(row_, stride_, fenv, _ts, _tr, _tkle)")
        return fn_text

    k0, k1 = find_function(kh, r"static KOKKOS_INLINE_FUNCTION auto kernel\(")
    extra = tab_fn + tabify(kh[k0:k1], "kernel") + "\n"
    pp = find_function(kh, r"static KOKKOS_INLINE_FUNCTION auto post\(")
    if pp:
        extra += tabify(kh[pp[0]:pp[1]], "post") + "\n"
    end = pp[1] if pp else k1
    new_kh = kh[:end] + "\n" + extra + kh[end:]

    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / "kernels.hh").write_text(new_khh)
    (a.out_dir / "kernel.hh").write_text(new_kh)


if __name__ == "__main__":
    main()
