#!/usr/bin/env python3
"""dagstat.py — static analyzer for NumTracer-emitted kernel headers.

Reconstructs the SSA DAG of every trace function in a *_kernels.hh (statements of the form
`const double sN = <op>;` as emitted by codegen/gen.hpp emit_rhs) and reports, per function
and in aggregate:

  - op counts by kind (const / var-load / add / sub / mul / neg) — the fp64 work proxy
  - max-live under the emitted order (register-pressure proxy)
  - fanout histogram, plus the two emission-opportunity counters:
      * fma-eligible: MULs whose single consumer is an ADD/SUB          (plan item B1)
      * single-use constants: RCONST slots referenced exactly once      (plan item B2)
  - return shape: real, complex two-sided, or complex half-zero ({0.0, sN} / {sN, 0.0})

With --support it additionally parses the companion *_kernel.hh (fill()/dr_/cosl1p/_interp
defs) and classifies every definition and every trace-fn SSA op by its transitive variable
support against the base kernel arguments {l1, cos1, cos2, phi, p, k}:

      K    : depends only on k                 — invariant across the whole integration
      PK   : depends only on {p, k}           — invariant per grid point (plan item G2)
      NODE : depends on a quadrature variable — genuinely per-node work

With --monomials it expands each trace function's return DAG to distinct monomials in the
f[i] env slots (capped; use for small/medium functions).

The ~146 KB `using _T = decltype(...)` line and any other pathological line are DROPPED
before parsing (see STALE_TEST_ORACLES.md — any grep over the raw file double-counts it).

Usage:
    dagstat.py ZA4_147_num_kernels.hh
    dagstat.py --support ZA3_gpu_kernel.hh ZA3_gpu_kernels.hh
    dagstat.py --monomials --funcs 'tr0|tr1' ZA3_147_num_kernels.hh
"""

import argparse
import collections
import re
import sys
from pathlib import Path

MAX_LINE = 100_000  # the _T decltype trap: drop, never parse
BASE_VARS = ("l1", "cos1", "cos2", "cos3", "phi", "p", "k")
NODE_VARS = {"l1", "cos1", "cos2", "cos3", "phi"}

STMT_RE = re.compile(r"^\s*const double (s\d+) = (.+);\s*$")
# The decorator run between `static` and the return type is whatever the emitter's
# "Decorator" option was set to, and it has grown: `static inline`, the CUDA
# `__host__ __device__` pair, the `__attribute__((noinline))` size gate, and the Kokkos
# macros used by every production (device) kernel. Match the run generically rather than
# enumerating orderings — an unmatched spelling makes this tool report ZERO functions and
# print "no trace functions parsed", which reads like a wrong-file mistake rather than a
# blind spot. Likewise `nt_complex_t` (the `NT_TRACE_COMPLEX` alias) superseded the literal
# `std::complex<double>` in the emitted signature, so complex traces went unseen everywhere.
DECOR_RE = (
    r"(?:__host__|__device__|__forceinline__|inline|constexpr"
    r"|__attribute__\(\(noinline\)\)|KOKKOS_[A-Z_]+)"
)
FUN_RE = re.compile(
    r"^static\s+(?:" + DECOR_RE + r"\s+)*"
    r"(void|double|float|nt_complex_t|std::complex<double>|auto)\s+(\w+)\("
)
# CrossTraceCSE emits ONE `void trace_all(const double *f, T *t)` per chunk instead of N `trN()`
# functions, and its results leave through `t[i] = ...;` stores rather than a `return`. Those stores
# are the roots — without them every op in the fused body looks dead and the reported cost is zero.
STORE_RE = re.compile(r"^\s*t\[\d+\] = (.+);\s*$")
RET_RE = re.compile(r"^\s*return\s+(.*);\s*$")
NUM_RE = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
IDENT_RE = re.compile(r"\b[A-Za-z_]\w*\b")
DEF_RE = re.compile(r"^\s*const (?:double|auto) (\w+) = (.+);\s*$")
FSLOT_RE = re.compile(r"^\s*f\[(\d+)\] = (.+);\s*$")


class Op:
    __slots__ = ("kind", "args", "value", "idx")

    def __init__(self, kind, args=(), value=None, idx=0):
        self.kind, self.args, self.value, self.idx = kind, args, value, idx


def read_filtered(path):
    dropped = 0
    lines = []
    for ln in Path(path).read_text(errors="replace").splitlines():
        if len(ln) > MAX_LINE:
            dropped += 1
            continue
        lines.append(ln)
    if dropped:
        print(f"[dagstat] {path}: dropped {dropped} pathological line(s) >{MAX_LINE} chars")
    return lines


# operand: a slot, an env load, or an inlined single-use constant `(lit)` (emitted parenthesised)
OPD = r"-?\s*(?:s\d+|f\[\d+\]|\((?:-?\d+\.?\d*|-?\.\d+)(?:[eE][-+]?\d+)?\))"
OPD_RE = re.compile(OPD)


def parse_rhs(rhs):
    """emit_stmt grammar: literal | f[N] | a*b | a+b | a-b | -a | fma(a,b,c) with operands
    sM, f[N], or a parenthesised inlined constant."""
    rhs = rhs.strip()
    if NUM_RE.match(rhs):
        return Op("const", value=float(rhs))
    if re.fullmatch(r"f\[\d+\]", rhs):
        return Op("var", args=(rhs,))

    def slots(*ops):
        return tuple(t for o in ops for t in re.findall(r"s\d+|f\[\d+\]", o))

    m = re.fullmatch(rf"({OPD})\s*([*+])\s*({OPD})", rhs)
    if m:
        return Op("mul" if m.group(2) == "*" else "add", args=slots(m.group(1), m.group(3)))
    m = re.fullmatch(rf"({OPD})\s*-\s*({OPD})", rhs)
    if m:
        return Op("sub", args=slots(m.group(1), m.group(2)))
    m = re.fullmatch(rf"-\s*({OPD})", rhs)
    if m:
        return Op("neg", args=slots(m.group(1)))
    m = re.fullmatch(rf"fma\(\s*({OPD})\s*,\s*({OPD})\s*,\s*({OPD})\s*\)", rhs)
    if m:
        return Op("fma", args=slots(m.group(1), m.group(2), m.group(3)))
    return Op("other", args=tuple(t for t in re.findall(r"s\d+|f\[\d+\]", rhs)))


def roots_of(fn):
    """Every result expression of a trace fn: its `return`, plus any fused `t[i] =` stores."""
    out = []
    if fn["ret"]:
        out.append(fn["ret"])
    out.extend(fn.get("stores", ()))
    return out


def parse_functions(lines):
    """{name: {"ops": {...}, "order": [...], "ret": <return str>, "stores": [<store str>...]}}."""
    fns = {}
    cur = None
    depth = 0
    for ln in lines:
        if cur is None:
            m = FUN_RE.match(ln)
            if m and m.group(2) not in ("powr", "fill"):
                cur = {"name": m.group(2), "ops": {}, "order": [], "ret": None, "stores": []}
                depth = ln.count("{") - ln.count("}")
                if depth <= 0:
                    depth = 1 if "{" in ln else 0
                    if depth == 0:
                        depth = 1  # opening brace on next line; close on matching
            continue
        depth += ln.count("{") - ln.count("}")
        m = STMT_RE.match(ln)
        if m:
            op = parse_rhs(m.group(2))
            op.idx = len(cur["order"])
            cur["ops"][m.group(1)] = op
            cur["order"].append(m.group(1))
            continue
        m = RET_RE.match(ln)
        if m:
            cur["ret"] = m.group(1)
        m = STORE_RE.match(ln)
        if m:
            cur["stores"].append(m.group(1))
        if depth <= 0:
            fns[cur["name"]] = cur
            cur = None
    return fns


def ret_slots(ret):
    return re.findall(r"s\d+", ret or "")


def ret_shape(ret):
    if ret is None:
        return "none"
    if re.fullmatch(r"s\d+", ret.strip()):
        return "real"
    z = re.search(r"\{\s*([^,{}]+)\s*,\s*([^,{}]+)\s*\}", ret)
    if z:
        a, b = z.group(1).strip(), z.group(2).strip()
        if NUM_RE.match(a) and float(a) == 0.0:
            return "imag-only"
        if NUM_RE.match(b) and float(b) == 0.0:
            return "real-only"
        return "two-sided"
    return "other"


def analyze_fn(fn):
    ops, order = fn["ops"], fn["order"]
    counts = collections.Counter(op.kind for op in ops.values())
    uses = collections.Counter()
    for op in ops.values():
        for a in op.args:
            uses[a] += 1
    for s in (x for r in roots_of(fn) for x in ret_slots(r)):
        uses[s] += 1
    # max-live under emitted order
    last_use = {}
    for slot in order:
        for a in ops[slot].args:
            if a.startswith("s"):
                last_use[a] = ops[slot].idx
    for s in (x for r in roots_of(fn) for x in ret_slots(r)):
        last_use[s] = len(order)
    live = maxlive = 0
    expiring = collections.defaultdict(list)
    for s, i in last_use.items():
        expiring[i].append(s)
    for slot in order:
        live += 1
        maxlive = max(maxlive, live)
        live -= len(expiring[ops[slot].idx])
    # B1: MULs with exactly one use, whose consumer is an ADD/SUB
    consumer = collections.defaultdict(list)
    for slot in order:
        for a in ops[slot].args:
            if a.startswith("s"):
                consumer[a].append(slot)
    fma_eligible = sum(
        1
        for slot in order
        if ops[slot].kind == "mul"
        and uses[slot] == 1
        and len(consumer[slot]) == 1
        and ops[consumer[slot][0]].kind in ("add", "sub")
    )
    # B2: constants used exactly once
    consts = [s for s in order if ops[s].kind == "const"]
    single_use_consts = sum(1 for s in consts if uses[s] == 1)
    arith = sum(counts[k] for k in ("add", "sub", "mul", "neg", "fma"))
    return {
        "counts": counts,
        "arith": arith,
        "maxlive": maxlive,
        "uses": uses,
        "fma_eligible": fma_eligible,
        "nconst": len(consts),
        "single_use_consts": single_use_consts,
        "distinct_consts": len({ops[s].value for s in consts}),
        "ret_shape": ret_shape(fn["ret"]),
    }


# ---------------------------------------------------------------- support (--support)

CALL_HEADS = re.compile(r"\b(\w+)\s*\(")
OPS_TOKENS = re.compile(r"[+\-*/]")


def expr_ops(expr):
    """Crude fp64-op weight of a free-form C++ expression: arithmetic tokens + calls."""
    return len(OPS_TOKENS.findall(expr)) + len(CALL_HEADS.findall(expr))


def build_support(kernel_lines, fslot_lines):
    """Support sets for named defs (cosl1p*, dr_*, _interp*, _cse*, ...) and f[] slots.

    Two passes: named defs first (dr_* feed the fill() call), then the f[] slots, then the
    named defs AGAIN so that anything referencing the filled env (`fenv`/`f`, e.g.
    `_interpN = ntRe(trK(fenv))`) inherits the union of the slot supports — a trace call is
    node-dependent even though `fenv` itself is not a base variable."""
    defs = {}      # name -> (support:set, ops:int, expr)
    order = []
    env_sup = set()

    def support_of(expr):
        sup = set()
        for t in IDENT_RE.findall(expr):
            if t in BASE_VARS:
                sup.add(t)
            elif t in defs:
                sup |= defs[t][0]
            elif t in ("fenv", "f"):
                sup |= env_sup
        return sup

    for ln in kernel_lines:
        m = DEF_RE.match(ln)
        if not m:
            continue
        name, expr = m.group(1), m.group(2)
        defs[name] = (support_of(expr), expr_ops(expr), expr)
        order.append(name)

    fslots = {}
    for ln in fslot_lines:
        m = FSLOT_RE.match(ln)
        if m:
            fslots[f"f[{m.group(1)}]"] = (support_of(m.group(2)), expr_ops(m.group(2)))
    env_sup = set().union(*(s for s, _ in fslots.values())) if fslots else set()
    for name in order:  # re-resolve with env_sup known
        defs[name] = (support_of(defs[name][2]), defs[name][1], defs[name][2])
    return defs, order, fslots


def classify(sup):
    if not (sup & NODE_VARS):
        if sup <= {"k"}:
            return "K"
        return "PK"
    return "NODE"


def support_report(kernel_lines, kernels_lines, fns):
    defs, order, fslots = build_support(kernel_lines, kernels_lines)
    wsum = collections.Counter()
    nsum = collections.Counter()
    k_calls = []
    for name in order:
        sup, w, expr = defs[name]
        cls = classify(sup)
        wsum[cls] += w
        nsum[cls] += 1
        if cls in ("K", "PK") and name.startswith("_interp"):
            k_calls.append((cls, name, expr if len(expr) < 90 else expr[:87] + "..."))
    tot = sum(wsum.values()) or 1
    print("\n--- support: named defs in kernel body (op-weighted) ---")
    for cls in ("K", "PK", "NODE"):
        print(f"  {cls:4} : {nsum[cls]:5d} defs   {wsum[cls]:6d} ops  ({100*wsum[cls]/tot:5.1f}%)")
    if k_calls:
        print(f"  hoistable _interp defs ({len(k_calls)}):")
        for cls, name, expr in k_calls:
            print(f"    [{cls:2}] {name} = {expr}")

    # trace-fn SSA ops by transitive support
    print("\n--- support: trace-fn SSA arithmetic by transitive support ---")
    gtot = collections.Counter()
    for fname, fn in sorted(fns.items()):
        slot_sup = {}
        cnt = collections.Counter()
        for slot in fn["order"]:
            op = fn["ops"][slot]
            if op.kind == "const":
                sup = frozenset()
            elif op.kind == "var":
                sup = frozenset(fslots.get(op.args[0], (set(), 0))[0])
            else:
                sup = frozenset().union(*(slot_sup.get(a) or frozenset(fslots.get(a, (set(), 0))[0]) for a in op.args)) if op.args else frozenset()
            slot_sup[slot] = sup
            if op.kind in ("add", "sub", "mul", "neg", "fma"):
                cnt[classify(set(sup))] += 1
        for cls in cnt:
            gtot[cls] += cnt[cls]
        tot = sum(cnt.values()) or 1
        if cnt["K"] + cnt["PK"]:
            print(f"  {fname:8} arith={tot:6d}  K={cnt['K']:5d} ({100*cnt['K']/tot:4.1f}%)  "
                  f"PK={cnt['PK']:5d} ({100*cnt['PK']/tot:4.1f}%)")
    tot = sum(gtot.values()) or 1
    print(f"  {'TOTAL':8} arith={tot:6d}  K={gtot['K']:5d} ({100*gtot['K']/tot:4.1f}%)  "
          f"PK={gtot['PK']:5d} ({100*gtot['PK']/tot:4.1f}%)  NODE={gtot['NODE']} ({100*gtot['NODE']/tot:4.1f}%)")


# ---------------------------------------------------------------- monomials (--monomials)

MONO_CAP = 5_000_000


def expand_monomials(fn):
    """Distinct monomials (in f[i]) of each returned slot; None if capped or non-DAG op hit."""
    memo = {}

    def poly(slot):
        if slot in memo:
            return memo[slot]
        if slot.startswith("f["):
            r = {(slot,): 1.0}
            memo[slot] = r
            return r
        op = fn["ops"][slot]
        if op.kind == "const":
            r = {(): op.value}
        elif op.kind == "var":
            r = {(op.args[0],): 1.0}
        elif op.kind == "neg":
            r = {m: -c for m, c in poly(op.args[0]).items()}
        elif op.kind in ("add", "sub"):
            a, b = poly(op.args[0]), poly(op.args[1])
            r = dict(a)
            sgn = 1.0 if op.kind == "add" else -1.0
            for m, c in b.items():
                r[m] = r.get(m, 0.0) + sgn * c
                if r[m] == 0.0:
                    del r[m]
        elif op.kind == "mul":
            a, b = poly(op.args[0]), poly(op.args[1])
            if len(a) * len(b) > MONO_CAP:
                raise OverflowError
            r = {}
            for ma, ca in a.items():
                for mb, cb in b.items():
                    key = tuple(sorted(ma + mb))
                    r[key] = r.get(key, 0.0) + ca * cb
                    if r[key] == 0.0:
                        del r[key]
        else:
            raise ValueError(f"non-polynomial op {op.kind} in {slot}")
        if len(r) > MONO_CAP:
            raise OverflowError
        memo[slot] = r
        return r

    try:
        polys = [poly(s) for r in roots_of(fn) for s in ret_slots(r)]
    except (OverflowError, ValueError) as e:
        return None, str(e)
    union = {}
    for p in polys:
        for m, c in p.items():
            union[m] = union.get(m, 0.0) + c
    return (sum(len(p) for p in polys), len(union)), None


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--support", action="store_true", help="variable-support classification (needs kernel.hh + kernels.hh)")
    ap.add_argument("--monomials", action="store_true")
    ap.add_argument("--funcs", default=None, help="regex filter on function names")
    ap.add_argument("--per-func", action="store_true", help="print one line per function")
    args = ap.parse_args()

    all_lines = []
    per_file = {}
    for f in args.files:
        per_file[f] = read_filtered(f)
        all_lines += per_file[f]

    fns = parse_functions(all_lines)
    if args.funcs:
        fns = {n: f for n, f in fns.items() if re.search(args.funcs, n)}
    if not fns:
        sys.exit("no trace functions parsed (wrong file? use the *_kernels.hh)")

    tot = collections.Counter()
    shapes = collections.Counter()
    tot_fma = tot_const = tot_single = tot_arith = 0
    rows = []
    for name, fn in sorted(fns.items(), key=lambda kv: (len(kv[0]), kv[0])):
        a = analyze_fn(fn)
        tot += a["counts"]
        shapes[a["ret_shape"]] += 1
        tot_fma += a["fma_eligible"]
        tot_const += a["nconst"]
        tot_single += a["single_use_consts"]
        tot_arith += a["arith"]
        rows.append((name, a))
        if args.per_func:
            c = a["counts"]
            print(f"{name:10} arith={a['arith']:6d} (mul={c['mul']:5d} add={c['add']:5d} sub={c['sub']:5d} "
                  f"neg={c['neg']:4d} fma={c['fma']:4d})  const={a['nconst']:5d} (1-use {a['single_use_consts']:5d}) "
                  f"maxlive={a['maxlive']:4d}  fma-elig={a['fma_eligible']:5d}  ret={a['ret_shape']}")

    print(f"\n=== {len(fns)} trace function(s), {' + '.join(Path(f).name for f in args.files)} ===")
    print(f"  arith ops        : {tot_arith:,}  (mul {tot['mul']:,} / add {tot['add']:,} / sub {tot['sub']:,} / neg {tot['neg']:,} / fma {tot['fma']:,})")
    print(f"  const slots      : {tot_const:,}  single-use {tot_single:,} ({100*tot_single/max(tot_const,1):.0f}%)   [B2 opportunity]")
    print(f"  fma-eligible MULs: {tot_fma:,} ({100*tot_fma/max(tot['mul'],1):.0f}% of muls)   [B1 opportunity]")
    print(f"  return shapes    : {dict(shapes)}")
    print(f"  max-live (worst) : {max((a['maxlive'] for _, a in rows), default=0)}")

    if args.monomials:
        print("\n--- monomial expansion ---")
        for name, fn in sorted(fns.items(), key=lambda kv: (len(kv[0]), kv[0])):
            r, err = expand_monomials(fn)
            if r is None:
                print(f"  {name:10} SKIPPED ({err or 'cap'})")
            else:
                print(f"  {name:10} monomials={r[0]:,} (distinct after union {r[1]:,})")

    if args.support:
        kernel_lines = []
        kernels_lines = []
        for f, lines in per_file.items():
            if f.endswith("_kernels.hh") or any(FSLOT_RE.match(l) for l in lines[:400]):
                kernels_lines += lines
            else:
                kernel_lines += lines
        support_report(kernel_lines, kernels_lines, fns)


if __name__ == "__main__":
    main()
