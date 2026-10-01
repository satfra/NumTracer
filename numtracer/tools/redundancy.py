#!/usr/bin/env python3
"""How redundant are the generator's sub-term contractions?

STALE (2026-08-08): this parser predates the Stage-4 pooled emission — dressed generators now
emit a chain POOL (`chp0()`) plus per-sub-term index arrays (`sdchR0()`/`sdslR0()`), so the
inline `dch<i>`/`dsl<i>` builders this script looks for no longer exist and it reports nothing
on regenerated sources. Its question (the dedup ceiling) was answered and the dedup LANDED
(see NUMTRACER_LEVER_B_BRIEF.md), so the parser was not ported. For kernel-side statistics use
tools/dagstat.py; if the redundancy question ever reopens, port load() to resolve the pools.

The generated main contracts one Dirac trace per (net i, sub-term j). The same trace recurs across
nets, so most of those contractions recompute something already computed. This counts how many
DISTINCT traces there are — the ceiling on what the sub-term dedup can save — and how often each is
reused, which is what the trace-cache policy is built on.

The key is (dnet, lnet, dch, dsl), matching what CodegenGenerator.m keys the dedup on. An earlier version of
this script keyed on (dch, dsl, lnet) only, omitting dnet; that undercounts distinct traces and so
overstates both the redundancy factor and the max reuse.

Usage:
    redundancy.py                      # every flow in numtracer/tests/gen
    redundancy.py <dir>                # every flow whose generator C++ is in <dir>
    redundancy.py <dir> <flow>         # one flow, e.g. gen_zaaqbq1_qcd_num
"""
import collections
import pathlib
import re
import statistics
import sys

HERE = pathlib.Path(__file__).resolve().parent
DEFAULT_GEN = HERE.parent / "tests" / "gen"


def split_top(s):
    """Split a braced-init element list on its top-level commas."""
    out, depth, cur = [], 0, []
    for ch in s:
        if ch in "{([":
            depth += 1
        elif ch in "})]":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    if "".join(cur).strip():
        out.append("".join(cur).strip())
    return out


def load(files):
    """Reassemble one flow's net builders -> {name: [element strings]}.

    Understands both forms ntChunkDef emits: the direct `RET name(){ return {E1, ...}; }`, and for
    oversized builders the `void name_cK(RET& o){ o.push_back(E); ... }` helpers, which have to be
    concatenated in chunk order to reproduce the element list.
    """
    src = "\n".join(f.read_text(errors="replace") for f in files)
    elems = {}
    for m in re.finditer(
        r"^\s*(?:std::vector<[^\n{]*?>)\s+(\w+)\(\)\s*\{\s*return\s*\{(.*?)\};\s*\}", src, re.S | re.M
    ):
        elems[m.group(1)] = split_top(m.group(2))

    chunks = collections.defaultdict(dict)
    for m in re.finditer(r"^\s*void\s+(\w+)_c(\d+)\([^)]*\)\s*\{(.*?)\}\s*$", src, re.S | re.M):
        chunks[m.group(1)][int(m.group(2))] = [
            e.strip() for e in re.findall(r"o\.push_back\((.*?)\);", m.group(3), re.S)
        ]
    for name, ck in chunks.items():
        elems[name] = [e for k in sorted(ck) for e in ck[k]]
    return elems


def traces_of(elems):
    """The (dn, ln, dch, dsl) key of every (net, sub-term), plus the per-net sub-term counts."""
    keys, per_net, i = [], [], 0
    while True:
        dn, ln = elems.get(f"dnet{i}"), elems.get(f"lnet{i}")
        if dn is None and ln is None:
            break
        dch, dsl = elems.get(f"dch{i}"), elems.get(f"dsl{i}")  # dressed path only
        n = len(ln) if ln is not None else len(dn)
        per_net.append(n)
        for j in range(n):
            keys.append((dn[j] if dn else "", ln[j] if ln else "",
                         dch[j] if dch else "", dsl[j] if dsl else ""))
        i += 1
    return keys, per_net


def report(flow, keys, per_net):
    if not keys:
        return
    c = collections.Counter(keys)
    tot, dist = len(keys), len(c)
    singles = sum(1 for v in c.values() if v == 1)

    print(f"\n=== {flow} ===")
    print(f"  contractions (net, sub-term) : {tot:,}")
    print(f"  distinct traces              : {dist:,}   redundancy {tot / dist:.2f}x")
    print(f"  refCount == 1 (singletons)   : {singles:,} ({100 * singles / dist:.1f}%)  -- caching one buys nothing")
    print(f"  refCount >= 2 (reused)       : {dist - singles:,} ({100 * (dist - singles) / dist:.1f}%)  -- the trace cache")
    print(f"  reuse histogram              : {dict(sorted(collections.Counter(c.values()).items())[:10])}")
    print(f"  nets                         : {len(per_net):,}")
    print(f"  sub-terms per net            : min={min(per_net)} "
          f"median={statistics.median(per_net):.0f} max={max(per_net)}")
    # Per-net scheduling is capped by the heaviest net: it alone cannot be split across threads.
    print(f"  heaviest net / mean net      : {max(per_net) / (tot / len(per_net)):.1f}x"
          f"  -- the load-balance argument")


def main():
    gen = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_GEN
    want = sys.argv[2] if len(sys.argv) > 2 else None
    if not gen.is_dir():
        sys.exit(f"no such directory: {gen}")

    flows = collections.defaultdict(list)
    for f in gen.glob("gen_*.cpp"):
        flows[re.sub(r"(_u\d+)?\.cpp$", "", f.name)].append(f)
    if not flows:
        sys.exit(f"no generator C++ (gen_*.cpp) in {gen}")

    for flow in sorted(flows):
        if want and flow != want:
            continue
        keys, per_net = traces_of(load(sorted(flows[flow])))
        report(flow, keys, per_net)


if __name__ == "__main__":
    main()
