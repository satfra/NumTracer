#!/usr/bin/env python3
"""noinline.py — apply the device __noinline__ size gate to a vendored kernels.hh.

gen.hpp's eff_decor() already implements this: on a device target, a trace function with more than
NT_GEN_NOINLINE_MIN (default 500) SSA instructions is emitted KOKKOS_FUNCTION +
__attribute__((noinline)) instead of KOKKOS_INLINE_FUNCTION. It is dead in both production codes —
`grep -c noinline` is 0 across all 16 flows — because the kernels predate the ntDeviceEnvPrefix fix
in Codegen*.m, so NT_GEN_DEVICE never reached the generator.

This script reproduces the gate's output on a vendored header so its effect can be measured per
architecture without a Wolfram regeneration. The threshold is a parameter precisely because 500 was
tuned on sm_89, where the 1:64 fp64 pipe hides spill; on a 1:2 part it is an open question.

Usage:
  noinline.py <in_kernels.hh> -o <out_kernels.hh> [--min 500]
  noinline.py <in_kernels.hh> --report
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

RE_FN = re.compile(r"^(\s*static\s+)(KOKKOS_INLINE_FUNCTION|KOKKOS_FORCEINLINE_FUNCTION|KOKKOS_FUNCTION)"
                   r"(\s+)(double|nt_complex_t)(\s+)(tr\d+)(\s*\()")
RE_SSA = re.compile(r"^\s*(?:\[\[maybe_unused\]\]\s*)?const\s+double\s+s\d+\s*=")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("header", type=Path)
    ap.add_argument("-o", "--out", type=Path)
    ap.add_argument("--min", type=int, default=500, help="instruction threshold (NT_GEN_NOINLINE_MIN)")
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()

    lines = args.header.read_text(errors="replace").splitlines()

    # First pass: measure every trace function.
    sizes: dict[int, tuple[str, int]] = {}
    i = 0
    while i < len(lines):
        m = RE_FN.match(lines[i])
        if not m:
            i += 1
            continue
        depth, seen, n, j = 0, False, 0, i
        while j < len(lines):
            depth += lines[j].count("{") - lines[j].count("}")
            if "{" in lines[j]:
                seen = True
            if RE_SSA.match(lines[j]):
                n += 1
            j += 1
            if seen and depth <= 0:
                break
        sizes[i] = (m.group(6), n)
        i = j

    # Second pass: rewrite the decorator of the ones over threshold.
    # KOKKOS_INLINE_FUNCTION expands to __forceinline__, which would contradict noinline, so the
    # gate swaps it for the plain KOKKOS_FUNCTION first — same order as eff_decor (gen.hpp:411-497).
    hit = 0
    out = list(lines)
    for idx, (name, n) in sizes.items():
        if n <= args.min:
            continue
        m = RE_FN.match(lines[idx])
        assert m
        out[idx] = (f"{m.group(1)}KOKKOS_FUNCTION __attribute__((noinline)){m.group(3)}"
                    f"{m.group(4)}{m.group(5)}{m.group(6)}{m.group(7)}"
                    + lines[idx][m.end():])  # keep the parameter list and everything after it
        hit += 1

    tot = sum(n for _, n in sizes.values())
    gated = sum(n for _, n in sizes.values() if n > args.min)
    print(f"{args.header.parent.name:<12} fns={len(sizes):<4} min={args.min:<5} "
          f"out-of-lined={hit:<4} ops_gated={gated}/{tot} ({100*gated//max(1,tot)}%)")

    if args.out and not args.report:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
