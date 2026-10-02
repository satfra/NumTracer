#!/usr/bin/env python3
"""attrprobe.py — attribute a flow's full-TU map-kernel fp64 SASS to buckets by cumulative stubbing.

The mkprobe.py census isolates the NumTracer emission (fill + traces); this tool measures the
COMPLEMENT, which FINDINGS.md's per-flow tables left as an unattributed 50-68% bucket. It compiles
the flow's real CT_map_1.cc four times, each time with kernel.hh progressively stubbed:

    full        the production kernel.hh, untouched (sanity: must reproduce the production count)
    nont        NumTracer traces stubbed (trN calls -> constants, fill dropped)
    nont_ni     + every interpolator lookup stubbed (index()/at()/direct calls -> constants)
    nont_ni_nr  + every regulator call (RB/RBdot/RF/RFdot) stubbed

Successive differences then give: NT-in-kernel, interpolator cost, regulator cost, and the final
残り = wrapper + angle algebra + weights + the kernel's combination polynomial.

The compile command is taken VERBATIM from the app build's compile_commands.json (same nvcc
wrapper, same flags, same -arch), with only the source path and -o redirected — so the numbers are
the production numbers. Counts are DADD|DMUL|DFMA|MUFU per SASS function; the map kernel is the
largest function in the TU (Kokkos emits it twice, local/global launch variants — take the max).

Usage: attrprobe.py <flow_dir> <build_dir_with_compile_commands> [-o results.csv] [--work DIR]
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

FP64 = re.compile(r"\b(DADD|DMUL|DFMA|DSETP|MUFU)\b")
FUNC = re.compile(r"Function : (\S+)")


def stub_nont(text: str) -> str:
    text = re.sub(r"[\w:]*::fill\(fenv[^;]*\);", "(void)fenv;", text)
    text = re.sub(r"[\w:]*::tr(\d+)\(fenv\)", lambda m: f"(1e-3 * ({m.group(1)} + 1))", text)
    # complex flows read the traces through ntRe/ntIm
    text = re.sub(r"const auto (_interp\d+) = nt(Re|Im)\([^;]*\);",
                  lambda m: f"const auto {m.group(1)} = 1e-3;", text)
    return text


def stub_nointerp(text: str) -> str:
    text = re.sub(r"const auto (_ix\d+) = \w+\.index\([^;]*\);",
                  lambda m: f"const auto {m.group(1)} = 30.5;", text)
    n = [0]

    def cst(m):
        n[0] += 1
        return f"const auto {m.group(1)} = {1.0 + 0.001 * n[0]:.6f};"

    text = re.sub(r"const auto (_interp\d+) = \w+\.at\(_ix\d+\);", cst, text)
    # direct dressing calls Name(arg) — everything that is NOT a regulator
    text = re.sub(r"const auto (_interp\d+) = (?!R[BF])\w+\((?!fenv)[^;]*\);", cst, text)
    return text


def stub_noreg(text: str) -> str:
    n = [0]

    def cst(m):
        n[0] += 1
        return f"const auto {m.group(1)} = {0.5 + 0.0001 * n[0]:.6f};"

    return re.sub(r"const auto (_interp\d+) = R[BF]\w*\([^;]*\);", cst, text)


VARIANTS = {
    "full": [],
    "nont": [stub_nont],
    "nont_ni": [stub_nont, stub_nointerp],
    "nont_ni_nr": [stub_nont, stub_nointerp, stub_noreg],
}


def fp64_of_biggest_fn(obj: Path) -> tuple[int, int]:
    sass = subprocess.run(["cuobjdump", "-sass", str(obj)], capture_output=True, text=True).stdout
    counts: dict[str, int] = {}
    fn = "?"
    for line in sass.splitlines():
        m = FUNC.search(line)
        if m:
            fn = m.group(1)
            counts.setdefault(fn, 0)
        elif FP64.search(line):
            counts[fn] = counts.get(fn, 0) + 1
    if not counts:
        return 0, 0
    top = sorted(counts.values(), reverse=True)
    return top[0], len(counts)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("flow_dir", type=Path)
    ap.add_argument("build_dir", type=Path)
    ap.add_argument("-o", "--out", type=Path)
    ap.add_argument("--work", type=Path, default=None)
    args = ap.parse_args()

    flow = args.flow_dir.resolve()
    db = json.load(open(args.build_dir / "compile_commands.json"))
    entry = next((e for e in db if str(Path(e["file"]).resolve()).endswith(f"{flow.name}/src/CT_map_1.cc")), None)
    if entry is None:
        print(f"no compile_commands entry for {flow.name}/src/CT_map_1.cc", file=sys.stderr)
        return 1

    work = (args.work or (Path(__file__).resolve().parent.parent / "variants" / "attr" / flow.name))
    rows = []
    for variant, stubs in VARIANTS.items():
        vdir = work / variant
        if vdir.exists():
            shutil.rmtree(vdir)
        shutil.copytree(flow, vdir)
        khh = vdir / "kernel.hh"
        text = khh.read_text()
        for s in stubs:
            text = s(text)
        khh.write_text(text)

        obj = vdir / "CT_map_1.o"
        src_orig = entry["file"]
        cmd = entry["command"]
        # compile_commands.json records the underlying host compiler; the build actually goes
        # through kokkos_launch_compiler + nvcc_wrapper (see the generated build.make rule).
        launcher = Path.home() / ".local/share/DiFfRG/bundled/bin/kokkos_launch_compiler"
        wrapper = Path.home() / ".local/share/DiFfRG/bundled/bin/nvcc_wrapper"
        if cmd.startswith("/usr/bin/c++ ") and launcher.exists():
            cmd = f"{launcher} {wrapper} /usr/bin/c++ /usr/bin/c++ " + cmd[len("/usr/bin/c++ "):]
        # redirect source and object; the original -o path stays valid but harmless to replace
        cmd = cmd.replace(src_orig, str(vdir / "src" / "CT_map_1.cc"))
        cmd = re.sub(r"-o \S+", f"-o {obj}", cmd)
        r = subprocess.run(cmd, shell=True, cwd=entry["directory"], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"{variant}: COMPILE FAILED\n{r.stderr[-2000:]}", file=sys.stderr)
            return 1
        fp64, nfn = fp64_of_biggest_fn(obj)
        rows.append((variant, fp64, nfn))
        print(f"{flow.name} {variant:>11}: map-kernel fp64 = {fp64}  ({nfn} SASS fns)")

    print(f"\n{flow.name} attribution (successive differences):")
    labels = ["NT in-kernel (traces/fill/combination)", "interpolator lookups", "regulators",
              "remainder (wrapper/angles/weights)"]
    prev = rows[0][1]
    for (variant, fp64, _), label in zip(rows[1:], labels):
        print(f"  {label:<42} {prev - fp64:>6}")
        prev = fp64
    print(f"  {labels[-1]:<42} {rows[-1][1]:>6}")

    if args.out:
        with open(args.out, "a") as f:
            for variant, fp64, nfn in rows:
                f.write(f"{flow.name},{variant},{fp64},{nfn}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
