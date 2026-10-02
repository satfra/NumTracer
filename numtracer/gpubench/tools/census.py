#!/usr/bin/env python3
"""census.py — L1 static census of a production flow (no compiler, no GPU).

Answers, per flow, the questions that decide which lever is worth building:

  traces   how many fp64 SSA ops the integrand costs per quadrature point, and how they are
           distributed over trace functions (the per-function histogram is what the device
           __noinline__ size gate keys on, so it sizes that lever directly)
  fill     how many transcendentals fill() issues, and how many DISTINCT ones it actually needs.
           The gap is the fill()-CSE lever: fill() is printed textually by mpoly_to_cpp and never
           goes through the hash-consed SSA the trace bodies use.
  wrapper  how many spline interpolator lookups the kernel class performs per quadrature point.
           Each carries a fp64 log1p through LogarithmicCoordinates1D::backward plus a global
           read, which is hidden under arithmetic on a 1:64-fp64 part and is not on a 1:2 part.

Usage:
  census.py <flow_dir> [<flow_dir> ...]        one directory per flow (contains kernels.hh, kernel.hh)
  census.py --csv results/l1.csv <flow_dir>...
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

# ---------------------------------------------------------------------------------------------
# The emitted forms we parse. These track gen.hpp's printer:
#   trace fn   "<decor> double trN([[maybe_unused]] const double *f) {"
#   SSA slot   "  const double sK = <expr>;"   (or "[[maybe_unused]] const double sK")
# ---------------------------------------------------------------------------------------------
RE_FN = re.compile(r"^\s*(?:static\s+)?\S.*?\b(?:double|nt_complex_t|void)\s+(tr\d+|fill)\s*\(")
RE_SSA = re.compile(r"^\s*(?:\[\[maybe_unused\]\]\s*)?const\s+double\s+s(\d+)\s*=\s*(.*?);\s*$")
RE_FMA = re.compile(r"\bfma\s*\(")
RE_CALL = re.compile(r"\b(sqrt|cos|sin|exp|log|log1p|pow|tanh|atan|cosh|sinh)\s*\(")
RE_POWR = re.compile(r"\bpowr\s*<\s*(-?\d+)\s*>")


def _brace_body(lines: list[str], start: int) -> tuple[list[str], int]:
    """Return the lines of the function body starting at `start` (the signature line)."""
    depth = 0
    body: list[str] = []
    i = start
    seen_open = False
    while i < len(lines):
        ln = lines[i]
        depth += ln.count("{") - ln.count("}")
        if "{" in ln:
            seen_open = True
        body.append(ln)
        i += 1
        if seen_open and depth <= 0:
            break
    return body, i


def _transcendental_census(body: list[str]) -> tuple[Counter, Counter]:
    """(total calls by name, distinct call-argument spellings by name).

    The distinct count is the honest lower bound on how many of these a perfect CSE would keep.
    We key on the *balanced* argument text so `sqrt(1. - powr<2>(cos1))` occurring 28 times counts
    once. Note this is a source-level count: nvcc's own GVN may already collapse identical calls,
    which is exactly why L2 (SASS) has to confirm before we change the emitter.
    """
    total: Counter = Counter()
    distinct: dict[str, set[str]] = {}
    text = "\n".join(body)
    for m in RE_CALL.finditer(text):
        name = m.group(1)
        # walk the balanced parenthesis to capture the argument text
        i = m.end() - 1
        depth = 0
        j = i
        while j < len(text):
            if text[j] == "(":
                depth += 1
            elif text[j] == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        arg = text[i + 1 : j]
        arg = re.sub(r"\s+", " ", arg).strip()
        total[name] += 1
        distinct.setdefault(name, set()).add(arg)
    return total, Counter({k: len(v) for k, v in distinct.items()})


def census_flow(flow_dir: Path) -> dict:
    kernels = flow_dir / "kernels.hh"
    wrapper = flow_dir / "kernel.hh"
    if not kernels.exists():
        raise SystemExit(f"{flow_dir}: no kernels.hh")

    lines = kernels.read_text(errors="replace").splitlines()

    fn_sizes: dict[str, int] = {}
    fill_body: list[str] = []
    i = 0
    while i < len(lines):
        m = RE_FN.match(lines[i])
        if not m:
            i += 1
            continue
        name = m.group(1)
        body, nxt = _brace_body(lines, i)
        if name == "fill":
            fill_body = body
        else:
            # a one-line forwarder (`trN(f){ return trK(f); }`) is a dedup alias, not real work
            n = sum(1 for ln in body if RE_SSA.match(ln))
            fn_sizes[name] = n
        i = nxt

    sizes = sorted(fn_sizes.values())
    total_ops = sum(sizes)
    fma = sum(len(RE_FMA.findall(ln)) for ln in lines)

    f_total, f_distinct = _transcendental_census(fill_body)
    fill_div = sum(ln.count("1.0/") for ln in fill_body)
    fill_powr = Counter(int(x) for ln in fill_body for x in RE_POWR.findall(ln))

    row = {
        "flow": flow_dir.name,
        "gen_lines": len(lines),
        "n_trace_fns": len(fn_sizes),
        "trace_ops": total_ops,
        "fn_median": sizes[len(sizes) // 2] if sizes else 0,
        "fn_max": sizes[-1] if sizes else 0,
        "fn_above_500": sum(1 for s in sizes if s > 500),
        "fma": fma,
        "nenv": _nenv(lines),
        "fill_div": fill_div,
        "fill_trans_total": sum(f_total.values()),
        "fill_trans_distinct": sum(f_distinct.values()),
        "fill_trans_waste": sum(f_total.values()) - sum(f_distinct.values()),
        "fill_powr_negative": sum(v for k, v in fill_powr.items() if k < 0),
    }
    row["_fill_detail"] = ", ".join(
        f"{k}:{f_total[k]}/{f_distinct[k]}" for k in sorted(f_total)
    )

    if wrapper.exists():
        wtext = wrapper.read_text(errors="replace")
        wlines = wtext.splitlines()
        row["interp_calls"] = len(set(re.findall(r"\b_interp\d+\b", wtext)))
        wt_total, wt_distinct = _transcendental_census(wlines)
        row["wrap_trans_total"] = sum(wt_total.values())
        row["wrap_trans_distinct"] = sum(wt_distinct.values())
        row["wrap_lines"] = len(wlines)
    else:
        row["interp_calls"] = row["wrap_trans_total"] = row["wrap_trans_distinct"] = 0
        row["wrap_lines"] = 0

    return row


def _nenv(lines: list[str]) -> int:
    for ln in lines:
        m = re.search(r"constexpr\s+int\s+nenv\s*=\s*(\d+)", ln)
        if m:
            return int(m.group(1))
    return 0


COLS = [
    "code", "flow", "gen_lines", "n_trace_fns", "trace_ops", "fn_median", "fn_max",
    "fn_above_500", "fma", "nenv", "interp_calls", "fill_div", "fill_trans_total",
    "fill_trans_distinct", "fill_trans_waste", "wrap_trans_total", "wrap_trans_distinct",
    "wrap_lines", "fill_powr_negative",
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--csv", type=Path, help="also write a csv here")
    ap.add_argument("--code", default="", help="label for the code the flows belong to (ym/nf2)")
    ap.add_argument("--detail", action="store_true", help="print the per-call fill breakdown")
    args = ap.parse_args()

    rows = []
    for d in args.dirs:
        if not (d / "kernels.hh").exists():
            continue
        r = census_flow(d)
        r["code"] = args.code or d.parent.parent.name
        rows.append(r)

    rows.sort(key=lambda r: -r["trace_ops"])

    hdr = f"{'code':>5} {'flow':<10} {'ops':>7} {'fns':>4} {'med':>5} {'max':>6} {'>500':>4} " \
          f"{'nenv':>4} {'interp':>6} {'div':>4} {'fillT':>6} {'fillD':>6} {'waste':>6} {'wrapT':>6}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['code']:>5} {r['flow']:<10} {r['trace_ops']:>7} {r['n_trace_fns']:>4} "
              f"{r['fn_median']:>5} {r['fn_max']:>6} {r['fn_above_500']:>4} {r['nenv']:>4} "
              f"{r['interp_calls']:>6} {r['fill_div']:>4} {r['fill_trans_total']:>6} "
              f"{r['fill_trans_distinct']:>6} {r['fill_trans_waste']:>6} {r['wrap_trans_total']:>6}")
        if args.detail and r["_fill_detail"]:
            print(f"{'':>5} {'':<10} fill total/distinct -> {r['_fill_detail']}")

    tot = sum(r["trace_ops"] for r in rows)
    print("-" * len(hdr))
    print(f"{'':>5} {'TOTAL':<10} {tot:>7}")

    neg = [r for r in rows if r["fill_powr_negative"]]
    if neg:
        print(f"\n!! negative powr<N> in fill(): {[r['flow'] for r in neg]} "
              f"— the emitted powr has no N<0 branch and returns 1.0 (the powr emission in CodegenGenerator.m)")

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"\ncsv -> {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
