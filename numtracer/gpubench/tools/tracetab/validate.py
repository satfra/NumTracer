#!/usr/bin/env python3
"""Stage 0c probe: is a k-coefficient table an accurate replacement for the inline trace?

For sampled nodes (qhat, angles, external point) the trace is expanded once in k (l1 = k*qhat),
with the inv slots either kept as runtime factors ("multi") or cleared into a univariate
numerator N(k) ("uni").  Coefficients are stored in float32 and the reconstruction is evaluated
in float32 at k over the whole flow range, against a float64 inline reference.  The current
float32 inline evaluation is the noise floor.

Error metric per (external point, k): max over loop nodes |T - T_ref| / max over loop nodes |T_ref|,
i.e. relative to the integrand scale, never pointwise.
"""
import math
import random
import re
import sys
from collections import defaultdict

import numpy as np

import analyse
from analyse import P

f32 = np.float32


def fma64(a, b, c):
    return a * b + c


def fma32(a, b, c):
    return f32(f32(a) * f32(b)) + f32(c)


def compile_trace(body, ret):
    lines = ["def tr(F, fma):"]
    for line in body.split("\n"):
        line = line.strip().replace("[[maybe_unused]] ", "")
        if not line.startswith("const "):
            continue
        lhs, rhs = line.split(" ", 2)[2].rstrip(";").split(" = ", 1)
        rhs = analyse.num.sub(r"\1", rhs)
        rhs = re.sub(r"f\[(\d+)\]", r"F[\1]", rhs)
        lines.append(f"  {lhs} = {rhs}")
    lines.append(f"  return {ret}")
    ns = {}
    exec(compile("\n".join(lines), "<tr>", "exec"), ns)
    return ns["tr"]


def compile_fill(src):
    m = re.search(r"void fill\(\w+ \*f, ([^)]*)\) \{(.*?)\n\}", src, re.S)
    args = [a.split()[-1] for a in m.group(1).split(",")]
    body = []
    for i, rhs in re.findall(r"f\[(\d+)\] = (.*?);\s*$", m.group(2), re.M):
        rhs = analyse.num.sub(r"\1", rhs)
        rhs = re.sub(r"powr<(-?\d+)>", r"pw\1", rhs).replace("pw-", "pwm")
        body.append(f"  f[{i}] = {rhs}")
    n = 1 + max(int(i) for i, _ in re.findall(r"f\[(\d+)\] = (.*?);\s*$", m.group(2), re.M))
    text = f"def fill({', '.join(args)}):\n  f = [0.0]*{n}\n" + "\n".join(body) + "\n  return f"
    ns = {"sqrt": math.sqrt, "cos": math.cos, "sin": math.sin}
    for e in range(-4, 9):
        ns[f"pw{e}" if e >= 0 else f"pwm{-e}"] = (lambda e: lambda x: x ** e)(e)
    exec(compile(text, "<fill>", "exec"), ns)
    return args, ns["fill"]


def main(path, nnode_ext=24, nnode_loop=24, nk=13, seed=3):
    src = open(path).read()
    slots, traces = analyse.parse(path)
    rad, inv, other = analyse.classify(slots)
    args, fill = compile_fill(src)
    trs = [(n, compile_trace(b, r)) for n, b, r in traces]
    P.nv = 1 + len(inv)
    sidx = {j: 1 + n for n, j in enumerate(inv)}
    rnd = random.Random(seed)
    ks = np.geomspace(0.034, 1000.0, nk)

    def ext_point():
        d = {}
        for a in args:
            if a in ("p", "S0"):
                d[a] = math.exp(rnd.uniform(math.log(0.005), math.log(250.0)))
            elif a == "S1":
                d[a] = rnd.uniform(0.0, 0.99)
            elif a in ("SPhi", "phi"):
                d[a] = rnd.uniform(0.0, 2 * math.pi)
        return d

    def loop_point():
        return {"qhat": rnd.uniform(0.01, 2.5), "cos1": rnd.uniform(-0.99, 0.99),
                "cos2": rnd.uniform(-0.99, 0.99), "phi": rnd.uniform(0, 2 * math.pi)}

    worst = defaultdict(lambda: {"inline32": 0.0, "multi": 0.0, "uni": 0.0})
    for _ in range(nnode_ext):
        ext = ext_point()
        loops = [loop_point() for _ in range(nnode_loop)]
        # per (trace, k): lists over loop nodes
        acc = defaultdict(lambda: defaultdict(list))
        for lp in loops:
            node = dict(ext)
            for a in ("cos1", "cos2", "phi"):
                if a in args and a not in ext:
                    node[a] = lp[a]

            def slots_at(k):
                kw = dict(node)
                kw["l1"] = k * lp["qhat"]
                return fill(**{a: kw[a] for a in args})

            # symbolic env: l1 -> qhat * k, inv symbolic, rest numeric (from fill at any k)
            base = slots_at(1.0)
            env = {}
            for j in range(len(base)):
                if j == rad[0]:
                    env[j] = P.var(0, lp["qhat"])
                elif j in sidx:
                    env[j] = P.var(sidx[j])
                else:
                    env[j] = base[j]
            # Q_j(k) = 1/inv_j(k) is a quadratic in k: recover it from 3 samples
            Q = {}
            for j in inv:
                y = [1.0 / slots_at(kk)[j] for kk in (0.0, 1.0, 2.0)]
                c0 = y[0]; c2 = (y[2] - 2 * y[1] + y[0]) / 2; c1 = y[1] - c0 - c2
                Q[j] = [c0, c1, c2]
            for name, tr in trs:
                poly = tr(env, fma64)
                live = poly.live() if isinstance(poly, P) else {(0,) * P.nv: float(poly)}
                maxpow = {j: max(m[sidx[j]] for m in live) for j in inv}
                # univariate numerator coefficients
                N = defaultdict(float)
                for m, c in live.items():
                    terms = {m[0]: c}
                    for j in inv:
                        for _ in range(maxpow[j] - m[sidx[j]]):
                            nt = defaultdict(float)
                            for d, v in terms.items():
                                for e, q in enumerate(Q[j]):
                                    nt[d + e] += v * q
                            terms = nt
                    for d, v in terms.items():
                        N[d] += v
                deg = max(N) if N else 0
                Ncoef = [f32(N.get(d, 0.0)) for d in range(deg + 1)]
                mcoef = [(m, f32(c)) for m, c in live.items()]
                for ik, k in enumerate(ks):
                    s = slots_at(k)
                    ref = tr(s, fma64)
                    s32 = [f32(x) for x in s]
                    in32 = tr(s32, fma32)
                    k32 = f32(k)
                    # multivariate: sum c k^n prod inv^m, all float32
                    mv = f32(0.0)
                    for m, c in mcoef:
                        t = c * k32 ** m[0]
                        for j in inv:
                            if m[sidx[j]]:
                                t = t * s32[j] ** m[sidx[j]]
                        mv = mv + t
                    # univariate: Horner in k, times prod inv^maxpow
                    uv = f32(0.0)
                    for c in reversed(Ncoef):
                        uv = uv * k32 + c
                    for j in inv:
                        uv = uv * s32[j] ** maxpow[j]
                    acc[(name, ik)]["ref"].append(ref)
                    acc[(name, ik)]["inline32"].append(float(in32))
                    acc[(name, ik)]["multi"].append(float(mv))
                    acc[(name, ik)]["uni"].append(float(uv))
        for (name, ik), d in acc.items():
            ref = np.array(d["ref"])
            scale = np.max(np.abs(ref))
            if scale == 0:
                continue
            for meth in ("inline32", "multi", "uni"):
                e = np.max(np.abs(np.array(d[meth]) - ref)) / scale
                worst[name][meth] = max(worst[name][meth], e)
    flow = path.split("/")[-2]
    print(f"# {flow}: max over (ext, k) of max_loop|T-T_ref|/max_loop|T_ref|  (k in [0.034, 1000])")
    print(f"{'tr':6s} {'inline32':>10s} {'multi32':>10s} {'uni32':>10s}")
    for name, w in worst.items():
        print(f"{name:6s} {w['inline32']:10.2e} {w['multi']:10.2e} {w['uni']:10.2e}")


if __name__ == "__main__":
    main(sys.argv[1], *[int(x) for x in sys.argv[2:]])
