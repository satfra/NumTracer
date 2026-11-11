#!/usr/bin/env python3
"""Stage 0c, homogeneous form.

Every trace is homogeneous of degree d in (l1, S0).  With l1 = k*qhat and the inv slots cleared,
N(k, S0) = T * prod_j Q_j^{m_j} is homogeneous of degree D = d + 2*sum m_j, so
    N = sum_n c_n k^n S0^(D-n),   c_n = c_n(qhat, angles, shape)   -- no S0 dependence.
Reconstruction with s = max(k, S0), r = min/max in [0, 1]:
    T = s^d * prod_j (s^2 inv_j)^{m_j} * Horner_r(c)
Coefficients are computed in float64 (once, at S0 = 1), stored and evaluated in DTYPE.
"""
import math
import random
import sys
from collections import defaultdict

import numpy as np
from fractions import Fraction
EXACT = '--exact' in sys.argv


def fmaq(a, b, c):
    return a * b + c

import analyse
from analyse import P
from validate import compile_fill, compile_trace, fma64, fma32

DTYPE = np.float32 if "--f64" not in sys.argv else np.float64


def main(path, n_shape=12, n_loop=16, n_s0=10, n_k=13, seed=7):
    src = open(path).read()
    slots, traces = analyse.parse(path)
    rad, inv, other = analyse.classify(slots)
    args, fill = compile_fill(src)
    trs = [(n, compile_trace(b, r)) for n, b, r in traces]
    P.nv = 1 + len(inv)
    sidx = {j: 1 + n for n, j in enumerate(inv)}
    sname = "S0" if "S0" in args else "p"
    rnd = random.Random(seed)
    ks = np.geomspace(0.034, 1000.0, n_k)
    s0s = np.geomspace(0.005, 250.0, n_s0)
    worst = defaultdict(lambda: {"inline32": 0.0, "hom": 0.0, "inline64": 0.0})

    # homogeneity degree of each trace: mode over random generic points
    degree = {}
    for name, tr in trs:
        votes = defaultdict(int)
        for _ in range(12):
            kw = {"S1": rnd.uniform(0, .99), "SPhi": rnd.uniform(0, 6.28), "qhat": rnd.uniform(.1, 2.5),
                  "cos1": rnd.uniform(-.99, .99), "cos2": rnd.uniform(-.99, .99), "phi": rnd.uniform(0, 6.28)}
            a1 = {**kw, "l1": kw["qhat"], sname: 1.0}; a2 = {**kw, "l1": 2 * kw["qhat"], sname: 2.0}
            t1 = tr(fill(**{a: a1[a] for a in args}), fma64); t2 = tr(fill(**{a: a2[a] for a in args}), fma64)
            if t1 != 0 and t2 != 0:
                votes[round(math.log2(abs(t2 / t1)))] += 1
        degree[name] = max(votes, key=votes.get) if votes else 0

    # mass dimension of each Q_j = 1/inv_j (2, or 0 for shape-only denominators)
    kw = {"S1": .4, "SPhi": 1.1, "qhat": .7, "cos1": .3, "cos2": -.2, "phi": .5}
    fa = fill(**{a: {**kw, "l1": .7, sname: 1.0}[a] for a in args})
    fb = fill(**{a: {**kw, "l1": 1.4, sname: 2.0}[a] for a in args})
    w = {j: round(math.log2(fa[j] / fb[j])) for j in inv}

    for _ in range(n_shape):
        shape = {"S1": rnd.uniform(0.0, 0.99), "SPhi": rnd.uniform(0, 2 * math.pi)}
        loops = [{"qhat": rnd.uniform(0.01, 2.5), "cos1": rnd.uniform(-0.99, 0.99),
                  "cos2": rnd.uniform(-0.99, 0.99), "phi": rnd.uniform(0, 2 * math.pi)} for _ in range(n_loop)]
        acc = defaultdict(lambda: defaultdict(list))
        for lp in loops:
            def kwargs(k, s0):
                kw = {**shape, **lp, "l1": k * lp["qhat"], sname: s0}
                return {a: kw[a] for a in args}

            base = fill(**kwargs(1.0, 1.0))
            env = {}
            for j in range(len(base)):
                env[j] = P.var(0, lp["qhat"]) if j == rad[0] else P.var(sidx[j]) if j in sidx else base[j]
            Q = {}
            for j in inv:
                y = [1.0 / fill(**kwargs(kk, 1.0))[j] for kk in (0.0, 1.0, 2.0)]
                c0 = y[0]; c2 = (y[2] - 2 * y[1] + y[0]) / 2; c1 = y[1] - c0 - c2
                Q[j] = [c0, c1, c2]
            for name, tr in trs:
                poly = tr(env, fma64)
                live = poly.live() if isinstance(poly, P) else {(0,) * P.nv: float(poly)}
                m = {j: max(mm[sidx[j]] for mm in live) for j in inv}
                N = defaultdict(float)
                for mono, c in live.items():
                    terms = {mono[0]: c}
                    for j in inv:
                        for _ in range(m[j] - mono[sidx[j]]):
                            nt = defaultdict(float)
                            for dd, v in terms.items():
                                for e, q in enumerate(Q[j]):
                                    nt[dd + e] += v * q
                            terms = nt
                    for dd, v in terms.items():
                        N[dd] += v
                # trace degree d from homogeneity (scale everything by 2)
                d = degree[name]
                D = d + sum(w[j] * m[j] for j in inv)
                c = [DTYPE(N.get(n, 0.0)) for n in range(D + 1)]
                for s0 in s0s:
                    for ik, k in enumerate(ks):
                        sl = fill(**kwargs(k, s0))
                        ref = float(tr([Fraction(x) for x in sl], fmaq)) if EXACT else tr(sl, fma64)
                        acc[(name, s0, ik)]['inline64'].append(tr(sl, fma64))
                        sl32 = [np.float32(x) for x in sl]
                        in32 = tr(sl32, fma32)
                        k_, s0_ = DTYPE(k), DTYPE(s0)
                        s = max(k_, s0_)
                        r = min(k_, s0_) / s
                        h = DTYPE(0)
                        # k <= S0: sum c_n r^n (r = k/S0); k > S0: sum c_n r^(D-n) (r = S0/k)
                        seq = reversed(c) if k_ <= s0_ else c
                        for cn in seq:
                            h = h * r + cn
                        fac = s ** d
                        for j in inv:
                            fac = fac * (s ** w[j] * DTYPE(sl[j])) ** m[j]
                        acc[(name, s0, ik)]["ref"].append(ref)
                        acc[(name, s0, ik)]["inline32"].append(float(in32))
                        acc[(name, s0, ik)]["hom"].append(float(h * fac))
        for key, dct in acc.items():
            ref = np.array(dct["ref"])
            sc = np.max(np.abs(ref))
            if sc == 0 or not np.isfinite(sc):
                continue
            for meth in ("inline32", "hom", "inline64"):
                e = np.max(np.abs(np.array(dct[meth]) - ref)) / sc
                worst[key[0]][meth] = max(worst[key[0]][meth], e)
    flow = path.split("/")[-2]
    tag = "f32" if DTYPE is np.float32 else "f64"
    print(f"# {flow}: homogeneous-table reconstruction ({tag}) vs {"exact" if EXACT else "float64 inline"}; max over (shape, S0, k)")
    print(f"{'tr':6s} {'inline32':>10s} {'hom':>10s} {'inline64':>10s}")
    for name, w in worst.items():
        print(f"{name:6s} {w['inline32']:10.2e} {w['hom']:10.2e} {w['inline64']:10.2e}")
    print(f"MAX    {max(w['inline32'] for w in worst.values()):10.2e} {max(w['hom'] for w in worst.values()):10.2e} {max(w['inline64'] for w in worst.values()):10.2e}")


if __name__ == "__main__":
    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(pos[0], *[int(x) for x in pos[1:]])
