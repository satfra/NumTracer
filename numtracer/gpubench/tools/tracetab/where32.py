import math, random, sys
from collections import defaultdict
from fractions import Fraction
import numpy as np
import analyse
from validate import compile_fill, compile_trace, fma64, fma32
def fmaq(a, b, c): return a * b + c
path = sys.argv[1]; src = open(path).read(); slots, traces = analyse.parse(path)
args, fill = compile_fill(src); trs = [(n, compile_trace(b, r)) for n, b, r in traces]
rnd = random.Random(11); sname = "S0" if "S0" in args else "p"
err = defaultdict(float)
for _ in range(6):
    shape = {"S1": rnd.uniform(0, .99), "SPhi": rnd.uniform(0, 6.28)}
    loops = [{"qhat": rnd.uniform(.01, 2.5), "cos1": rnd.uniform(-.99, .99), "cos2": rnd.uniform(-.99, .99), "phi": rnd.uniform(0, 6.28)} for _ in range(6)]
    for s0 in np.geomspace(0.005, 250, 9):
        for k in np.geomspace(0.034, 1000, 9):
            for name, tr in trs:
                ref = []; i32 = []
                for lp in loops:
                    kw = {**shape, **lp, "l1": k * lp["qhat"], sname: s0}
                    sl = fill(**{a: kw[a] for a in args})
                    ref.append(float(tr([Fraction(x) for x in sl], fmaq)))
                    i32.append(float(tr([np.float32(x) for x in sl], fma32)))
                sc = max(abs(x) for x in ref)
                if sc == 0: continue
                e = max(abs(a - b) for a, b in zip(i32, ref)) / sc
                b = math.floor(math.log10(s0 / k))
                err[b] = max(err[b], e)
print(path.split('/')[-2], " ".join(f"S0/k~1e{b}:{err[b]:.0e}" for b in sorted(err)))
