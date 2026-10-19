#!/usr/bin/env python3
"""fillcse.py — source-rewrite prototype of the fill()-CSE lever.

fill() is the one part of the emitted kernel that never goes through NumTracer's hash-consed SSA:
trace bodies are lowered by real_cse.hpp's RBuilder, while fill() is printed textually by
mpoly_to_cpp (poly_cpp.hpp:46). The result is visible in the production kernels — nf2 ZA4
issues 51 transcendental calls for 4 distinct values, ZAAqbq1/2 issue 103 for 4 — and the same
subexpression that was just stored into f[4] is re-substituted verbatim inside f[1], f[2], f[8]...

This script rewrites a vendored kernels.hh into an equivalent one with fill() CSE'd, so the lever
can be measured (L2 SASS diff, L3 runtime) BEFORE any emitter change is written. Two passes:

  1. slot reuse   — if f[i]'s RHS contains the exact text of an already-assigned f[j], use f[j].
  2. hoisting     — repeatedly hoist the shortest balanced call expression that occurs more than
                    once into a `const double _t<k>`, substituting everywhere. Shortest-first makes
                    nesting fall out naturally: powr<2>(cos1) is hoisted before sqrt(1. - <that>).

Substitution is precedence-safe by construction: a complete call expression is always replaced by a
single identifier, and the emitter already parenthesises the slot expressions it inlines.

Usage:
  fillcse.py <in_kernels.hh> -o <out_kernels.hh>          rewrite
  fillcse.py <in_kernels.hh> --report                     just count what it would save
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

RE_FILL_SIG = re.compile(r"^.*\bvoid\s+fill\s*\(")
RE_ASSIGN = re.compile(r"^(\s*)f\[(\d+)\]\s*=\s*(.*?);\s*$")
RE_CALL_START = re.compile(r"\b(sqrt|cos|sin|exp|log|log1p|tanh|atan|cosh|sinh|powr\s*<\s*-?\d+\s*>)\s*\(")

MIN_SLOT_REUSE_LEN = 20  # never substitute a trivial RHS like `p` or `l1` as a substring


def _balanced(text: str, open_idx: int) -> int:
    """Index just past the ')' matching the '(' at open_idx."""
    depth = 0
    i = open_idx
    while i < len(text):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return -1


def _call_exprs(text: str) -> list[str]:
    """Every complete call expression in `text`, innermost included."""
    out = []
    for m in RE_CALL_START.finditer(text):
        op = text.index("(", m.end() - 1)
        end = _balanced(text, op)
        if end > 0:
            out.append(text[m.start():end])
    return out


def cse_fill(body: list[str]) -> tuple[list[str], dict]:
    """Rewrite the fill() body. Returns (new_body, stats)."""
    head: list[str] = []
    assigns: list[tuple[str, int, str]] = []  # (indent, slot, rhs)
    tail: list[str] = []

    seen_assign = False
    for ln in body:
        m = RE_ASSIGN.match(ln)
        if m:
            seen_assign = True
            assigns.append((m.group(1), int(m.group(2)), m.group(3)))
        elif not seen_assign:
            head.append(ln)
        else:
            tail.append(ln)

    before_calls = sum(len(_call_exprs(r)) for _, _, r in assigns)

    # ---- pass 1: reuse already-assigned slots -------------------------------------------------
    reuse = 0
    for i in range(len(assigns)):
        ind, slot, rhs = assigns[i]
        for j in range(i):
            _, pslot, prhs = assigns[j]
            if len(prhs) < MIN_SLOT_REUSE_LEN:
                continue
            if prhs in rhs:
                n = rhs.count(prhs)
                rhs = rhs.replace(prhs, f"f[{pslot}]")
                reuse += n
        assigns[i] = (ind, slot, rhs)

    # ---- pass 2: hoist repeated call expressions ----------------------------------------------
    temps: list[str] = []
    hoisted = 0
    while True:
        texts = temps + [r for _, _, r in assigns]
        counts: dict[str, int] = {}
        for t in texts:
            for e in _call_exprs(t):
                counts[e] = counts.get(e, 0) + 1
        cands = [e for e, c in counts.items() if c >= 2]
        if not cands:
            break
        # shortest first: hoisting the innermost repeated expression lets the outer ones collapse
        pick = min(cands, key=lambda e: (len(e), e))
        name = f"_t{len(temps)}"
        temps = [t.replace(pick, name) for t in temps]
        temps.append(pick)
        assigns = [(i, s, r.replace(pick, name)) for i, s, r in assigns]
        # the temp we just appended must keep its own definition intact
        temps[-1] = pick
        hoisted += 1
        if hoisted > 4096:  # runaway guard
            break

    new_body = list(head)
    for k, t in enumerate(temps):
        new_body.append(f"  const double _t{k} = {t};")
    for ind, slot, rhs in assigns:
        new_body.append(f"{ind}f[{slot}] = {rhs};")
    new_body.extend(tail)

    after_calls = sum(len(_call_exprs(t)) for t in temps) + \
        sum(len(_call_exprs(r)) for _, _, r in assigns)

    return new_body, {
        "slots": len(assigns),
        "calls_before": before_calls,
        "calls_after": after_calls,
        "slot_reuses": reuse,
        "temps": len(temps),
    }


def process(path: Path) -> tuple[str, dict]:
    lines = path.read_text(errors="replace").splitlines()
    start = None
    for i, ln in enumerate(lines):
        if RE_FILL_SIG.match(ln):
            start = i
            break
    if start is None:
        raise SystemExit(f"{path}: no fill() found")

    depth = 0
    seen = False
    end = start
    for i in range(start, len(lines)):
        depth += lines[i].count("{") - lines[i].count("}")
        if "{" in lines[i]:
            seen = True
        if seen and depth <= 0:
            end = i
            break

    body = lines[start:end + 1]
    new_body, stats = cse_fill(body)
    out = lines[:start] + new_body + lines[end + 1:]
    return "\n".join(out) + "\n", stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("header", type=Path)
    ap.add_argument("-o", "--out", type=Path)
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()

    text, st = process(args.header)
    print(f"{args.header.parent.name:<12} slots={st['slots']:<3} "
          f"calls {st['calls_before']:>4} -> {st['calls_after']:<4} "
          f"(slot-reuses {st['slot_reuses']}, temps {st['temps']})")
    if args.out and not args.report:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
