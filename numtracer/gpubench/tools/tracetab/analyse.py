#!/usr/bin/env python3
"""Stage 0b probe: how big would a k-polynomial trace table be?

For every trN in an emitted *_kernels.hh, substitute l1 = k*qhat, keep the inv slots symbolic,
put random numbers into every other slot (angles, external momenta), and expand.  The result is
T = sum_{n,m} C_{n,m} k^n prod inv^m with node-fixed C: the number of monomials E is the number
of table entries per node.  Also reports the univariate form (inv slots expressed through k,
denominators cleared): E_uni = deg_k(N) + 1.

Float literals make exact cancellations inexact, so each monomial carries its absolute mass and
counts as zero when |c| < TOL * mass.
"""
import random
import re
import sys
from collections import defaultdict

TOL = 1e-5


class P:
    """Sparse polynomial: exponent tuple -> [coef, absmass]."""
    __slots__ = ("t",)
    nv = 0

    def __init__(self, t=None):
        self.t = t if t is not None else {}

    @staticmethod
    def const(c):
        return P({(0,) * P.nv: [c, abs(c)]})

    @staticmethod
    def var(i, c=1.0):
        e = [0] * P.nv
        e[i] = 1
        return P({tuple(e): [c, abs(c)]})

    def __add__(self, o):
        o = o if isinstance(o, P) else P.const(o)
        t = {k: v[:] for k, v in self.t.items()}
        for k, (c, a) in o.t.items():
            if k in t:
                t[k][0] += c
                t[k][1] += a
            else:
                t[k] = [c, a]
        return P(t)

    __radd__ = __add__

    def __neg__(self):
        return P({k: [-c, a] for k, (c, a) in self.t.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, P) else P.const(o)))

    def __rsub__(self, o):
        return P.const(o) - self

    def __mul__(self, o):
        if not isinstance(o, P):
            return P({k: [c * o, a * abs(o)] for k, (c, a) in self.t.items()})
        t = {}
        for k1, (c1, a1) in self.t.items():
            for k2, (c2, a2) in o.t.items():
                k = tuple(x + y for x, y in zip(k1, k2))
                if k in t:
                    t[k][0] += c1 * c2
                    t[k][1] += a1 * a2
                else:
                    t[k] = [c1 * c2, a1 * a2]
        return P(t)

    __rmul__ = __mul__

    def live(self):
        return {k: c for k, (c, a) in self.t.items() if a > 0 and abs(c) > TOL * a}


def fma(a, b, c):
    return a * b + c


def parse(path):
    src = open(path).read()
    fill = re.search(r"void fill\(([^)]*)\) \{(.*?)\n\}", src, re.S)
    slots = {}
    for i, rhs in re.findall(r"f\[(\d+)\] = (.*?);\s*$", fill.group(2), re.M):
        slots[int(i)] = rhs.strip()
    traces = []
    for m in re.finditer(r"\b(tr\d+)\(\[\[maybe_unused\]\] const \w+ \*f\) \{\n(.*?)\n\}", src, re.S):
        body, ret = ("\n" + m.group(2)).rsplit("\n  return ", 1)
        traces.append((m.group(1), body, num.sub(r"\1", ret.rstrip(";"))))
    return slots, traces


def classify(slots):
    """-> (radial slot ids, inv slot ids, other slot ids)."""
    rad = [i for i, r in slots.items() if r == "l1"]
    inv = [i for i, r in slots.items() if r.startswith("1.f/(") or r.startswith("1./(") or r.startswith("1.0/(")]
    other = [i for i in slots if i not in rad and i not in inv]
    return rad, inv, other


num = re.compile(r"(?<![\w.])(\d+\.?\d*(?:e[+-]?\d+)?)f\b")


def run_trace(body, ret, env):
    loc = {"fma": fma}
    nops = 0
    nflop = 0
    for line in body.split("\n"):
        line = line.strip().replace("[[maybe_unused]] ", "")
        if not line.startswith("const "):
            continue
        lhs, rhs = line[len("const float "):].rstrip(";").split(" = ", 1) if line.startswith("const float") else \
            line[len("const double "):].rstrip(";").split(" = ", 1)
        rhs = num.sub(r"\1", rhs)
        rhs = re.sub(r"f\[(\d+)\]", r"F[\1]", rhs)
        if re.fullmatch(r"-?\(?-?[\d.e+-]+\)?", rhs) is None and not rhs.startswith("F["):
            nops += 1
            nflop += 2 if rhs.startswith("fma(") else 1
        loc[lhs] = eval(rhs, {"F": env, "fma": fma}, loc)
    return eval(ret, {}, loc), nops, nflop


def univariate_degree(poly, inv_ids, slot_index):
    """deg_k of N = T * prod_j D_j^{maxpow_j}, D_j quadratic in k (generic)."""
    best = 0
    maxpow = defaultdict(int)
    for mono in poly:
        for j in inv_ids:
            maxpow[j] = max(maxpow[j], mono[slot_index[j]])
    for mono in poly:
        d = mono[0]
        for j in inv_ids:
            d += 2 * (maxpow[j] - mono[slot_index[j]])
        best = max(best, d)
    return best, sum(maxpow.values())


def main():
    rnd = random.Random(1)
    rows = []
    for path in sys.argv[1:]:
        slots, traces = parse(path)
        rad, inv, other = classify(slots)
        if len(rad) != 1:
            print(f"# {path}: radial slot not found ({rad}); skipped", file=sys.stderr)
            continue
        # variables: 0 = k, 1.. = inv slots
        P.nv = 1 + len(inv)
        slot_index = {j: 1 + n for n, j in enumerate(inv)}
        env = {}
        env[rad[0]] = P.var(0, rnd.uniform(0.3, 2.0))      # l1 = k * qhat
        for j in inv:
            env[j] = P.var(slot_index[j])
        for j in other:
            env[j] = rnd.uniform(-1.0, 1.0) if "cos" in slots[j] or "sqrt" in slots[j] else rnd.uniform(0.3, 2.0)
        flow = path.split("/")[-2]
        for name, body, ret in traces:
            val, nops, nflop = run_trace(body, ret, env)
            live = val.live() if isinstance(val, P) else {(0,) * P.nv: val}
            kdeg = max((m[0] for m in live), default=0)
            ninvmono = len({m[1:] for m in live})
            duni, sumpow = univariate_degree(live, inv, slot_index)
            rows.append((flow, name, nops, nflop, len(live), kdeg, ninvmono, duni + 1))
    print(f"{'flow':14s} {'tr':5s} {'ops':>6s} {'flop':>6s} {'E':>5s} {'deg_k':>5s} {'#inv-mono':>9s} {'E_uni':>5s}")
    for r in rows:
        print(f"{r[0]:14s} {r[1]:5s} {r[2]:6d} {r[3]:6d} {r[4]:5d} {r[5]:5d} {r[6]:9d} {r[7]:5d}")
    tot = defaultdict(lambda: [0, 0, 0, 0, 0])
    for r in rows:
        t = tot[r[0]]
        t[0] += 1; t[1] += r[2]; t[2] += r[3]; t[3] += r[4]; t[4] += r[7]
    print("\n# per flow: ntr, ops, flop, sum E, sum E_uni")
    for f, t in tot.items():
        print(f"{f:14s} {t[0]:4d} {t[1]:7d} {t[2]:7d} {t[3]:6d} {t[4]:6d}")


if __name__ == "__main__":
    main()
