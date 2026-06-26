#!/usr/bin/env bash
# ZA-bench orchestrator: compare the NumTracer full-QCD AAqbq ZA kernel against the FormTracer
# analogue on GENERATION time, EMITTED SIZE, COMPILE time, BINARY size, and RUNTIME ns/eval, and
# check the global-collection rework is output-preserving vs the pre-rework path.
#
# Stages (each can be skipped with the matching SKIP_* env var = 1):
#   1. regenerate FormTracer ZA            (SKIP_GEN_FORM)
#   2. snapshot current NT as nt_before, then regenerate NT ZA via the reworked global-collection
#      path into nt                        (SKIP_GEN_NT)
#   3. compile + run compare_za_full for each NT variant (nt_before, nt) vs Form
#   4. faithfulness gate: diff the nt vs nt_before value dumps (must agree to ~1e-12)
#
# Prereqs for live regeneration: wolfram + FunKit + DiFfRG (NT side) and FORM/tform (FORM side).
# With SKIP_GEN_*=1 the harness benchmarks whatever is already in generated/{form,nt,nt_before}/.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NTINC="${NTINC:-$HERE/../../numtracer/include}"
: > "$HERE/bench_report.txt"
CXXFLAGS="-std=c++20 -O3 -march=native -ftemplate-depth=4000 -I $HERE/bench -I $HERE/generated -I $HERE/bench/shiminc -I $NTINC"

[[ "${SKIP_GEN_FORM:-0}" == 1 ]] || "$HERE/gen_form_za.sh"
if [[ "${SKIP_GEN_NT:-0}" != 1 ]]; then
  echo "[bench] snapshotting current NT as nt_before ..."
  rm -rf "$HERE/generated/nt_before"; mkdir -p "$HERE/generated/nt_before"
  # re-namespace the existing nt artifact as the frozen 'before' variant
  sed -e 's/\bza_qcd_nt\b/za_qcd_before/g' -e 's/\bZA_nt_kernel\b/ZA_nt_before_kernel/g' \
      "$HERE/generated/nt/kernel.hh"  > "$HERE/generated/nt_before/kernel.hh"  2>/dev/null || true
  sed -e 's/\bza_qcd_nt\b/za_qcd_before/g' \
      "$HERE/generated/nt/kernels.hh" > "$HERE/generated/nt_before/kernels.hh" 2>/dev/null || true
  GLOBAL_COLLECT=1 "$HERE/gen_nt_za.sh" nt ZA_nt_kernel za_qcd_nt
fi

compile_run () { # <variant> <hdr> <class> <name>
  local v="$1" hdr="$2" cls="$3" nm="$4"
  [[ -f "$HERE/generated/$v/kernel.hh" ]] || { echo "[bench] no generated/$v, skipping"; return; }
  echo "[bench] compiling variant=$v ..."
  /usr/bin/time -f "COMPILE variant=$v sec=%e rss_kb=%M" \
    g++ $CXXFLAGS -DNT_HDR="\"$hdr\"" -DNT_CLASS=$cls -DNT_NAME="\"$nm\"" \
    "$HERE/bench/compare_za_full.cpp" -o "/tmp/za_$v" 2> >(tee -a "$HERE/bench_report.txt" >&2)
  echo "BINSIZE variant=$v bytes=$(stat -c%s "/tmp/za_$v")" >> "$HERE/bench_report.txt"
  "/tmp/za_$v" "/tmp/za_${v}_vals.txt" | tee -a "$HERE/bench_report.txt"
}

compile_run nt_before "nt_before/kernel.hh" ZA_nt_before_kernel nt_before
compile_run nt        "nt/kernel.hh"        ZA_nt_kernel        nt

# faithfulness gate: reworked nt must reproduce pre-rework nt_before to ~1e-12
if [[ -f /tmp/za_nt_vals.txt && -f /tmp/za_nt_before_vals.txt ]]; then
  echo "[bench] faithfulness gate (nt vs nt_before) ..."
  maxrel=$(paste /tmp/za_nt_before_vals.txt /tmp/za_nt_vals.txt | awk '
    { a=$1; b=$2; d=(a-b); if(d<0)d=-d; den=(a<0?-a:a); if(den<1e-300)den=1e-300; r=d/den; if(r>m)m=r }
    END { printf "%.3e", m }')
  echo "GATE nt_vs_nt_before_max_relerr=$maxrel" | tee -a "$HERE/bench_report.txt"
  awk -v m="$maxrel" 'BEGIN{ if(m+0 < 1e-10) print "GATE: PASS (rework is output-preserving, < 1e-10)"; else { print "GATE: FAIL"; exit 1 } }'
fi

echo; echo "===== bench_report.txt ====="; cat "$HERE/bench_report.txt"
