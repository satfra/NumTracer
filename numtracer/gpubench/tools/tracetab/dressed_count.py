#!/usr/bin/env python3
"""dressed_count.py -- how big would trace tables be for flows whose traces carry dressing slots?

A dressed trace is sum_alpha D^alpha * T_alpha(kinematics), D = the dressing slots (dr_*, fixed at run time).
Each T_alpha is a kinematic polynomial and can be tabulated like a dressing-free trace: after clearing its inv
denominators, the entries per node are the length of its live k-power range (l1 = k * qhat).

For each trace: expand with l1 -> k, the inv and dressing slots symbolic and every other slot a random number; group the
live monomials by dressing exponent; per group clear the denominators and count the live k-range.

Usage: dressed_count.py <kernels.hh>...
"""
import random
import re
import sys
import time
from collections import defaultdict

import analyse
from analyse import P
from validate import compile_fill, compile_trace, fma64


def main():
    rnd = random.Random(3)
    print(f"{'flow':16s} {'traces':>6s} {'flop/node':>10s} {'alpha':>7s} {'alpha distinct':>14s} {'entries E':>9s} "
          f"{'flop/E':>7s} {'sec':>5s}")
    for path in sys.argv[1:]:
        t0 = time.time()
        src = open(path).read()
        slots, traces = analyse.parse(path)
        rad = [i for i, r in slots.items() if r == "l1"]
        inv = [i for i, r in slots.items() if re.match(r"1\.0?f?/\(", r)]
        dress = [i for i, r in slots.items() if re.fullmatch(r"dr_\d+", r)]
        if len(rad) != 1:
            print(f"{path}: no radial slot; skipped")
            continue
        args, fill = compile_fill(src)
        sym = inv + dress
        P.nv = 1 + len(sym)
        sidx = {j: 1 + n for n, j in enumerate(sym)}
        ninv = len(inv)
        # one random generic node; dressing args are irrelevant (their slots are symbolic)
        kw = {a: rnd.uniform(0.3, 0.9) for a in args}
        vals = fill(**{**kw, "l1": 1.0})
        q = [fill(**{**kw, "l1": x}) for x in (0.0, 1.0, 2.0)]
        Q = {}
        for j in inv:
            y = [1.0 / q[i][j] for i in range(3)]
            c2 = (y[2] - 2 * y[1] + y[0]) / 2
            Q[j] = [y[0], y[1] - y[0] - c2, c2]
        env = {j: (P.var(0, 1.0) if j == rad[0] else P.var(sidx[j]) if j in sidx else vals[j]) for j in range(len(vals))}
        flop = 0
        n_alpha = 0
        alphas = set()
        E = 0
        for name, body, ret in traces:
            # pure-imaginary traces return nt_complex_t{0.0, s}: the structure is that of re + im
            ret = re.sub(r"nt_complex_t\{([^,]+), ([^}]+)\}", r"(\1) + (\2)", ret)
            flop += sum(2 if l.split(" = ", 1)[1].startswith("fma(") else 1
                        for l in body.split("\n") if " = " in l and not re.search(r"= -?\(?-?[\d.e+-]+f?\)?;$", l)
                        and "= f[" not in l)
            poly = compile_trace(body, ret)(env, fma64)
            if not isinstance(poly, P):
                continue
            groups = defaultdict(dict)
            for mono, (c, a) in poly.t.items():
                if a > 0 and abs(c) > 1e-6 * a:
                    groups[mono[1 + ninv:]][mono[: 1 + ninv]] = c
            for alpha, terms in groups.items():
                n_alpha += 1
                alphas.add(alpha)
                m = {j: max(t[sidx[j]] for t in terms) for j in inv}
                N = defaultdict(float)
                for t, c in terms.items():
                    cur = {t[0]: c}
                    for j in inv:
                        for _ in range(m[j] - t[sidx[j]]):
                            nxt = defaultdict(float)
                            for dd, v in cur.items():
                                for e, qe in enumerate(Q[j]):
                                    nxt[dd + e] += v * qe
                            cur = nxt
                    for dd, v in cur.items():
                        N[dd] += v
                scale = max(abs(v) for v in N.values()) if N else 0
                live = [dd for dd, v in N.items() if abs(v) > 1e-9 * scale]
                if live:
                    E += max(live) - min(live) + 1
        dt = time.time() - t0
        flow = path.split("/")[-2]
        print(f"{flow:16s} {len(traces):6d} {flop:10d} {n_alpha:7d} {len(alphas):14d} {E:9d} {flop / max(E, 1):7.0f} {dt:5.0f}",
              flush=True)


if __name__ == "__main__":
    main()
