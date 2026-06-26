#!/usr/bin/env bash
# Live-regenerate the FormTracer full-QCD AAqbq ZA kernel and copy it into generated/form/.
# Reports the FormTracer generation time. By default the FORM trace cache is emptied first so the
# trace is genuinely recomputed (else FormTracer's AutoSaveRestore serves the cached result and the
# time is ~0 — see memory `form-tracecache-regen`); pass KEEP_CACHE=1 to time a cache hit instead.
#
# NOTE: QCD.m is monolithic (it regenerates all 11 flows). We empty ONLY the ZA-relevant trace
# cache (FlowAA*) so re-tracing is scoped to ZA; the other flows are served from cache. The ZA
# MakeKernel time is also printed by FormTracer to the console log captured in form_gen.log.
#
# Usage:  ./gen_form_za.sh
# Env overrides:
#   QCD_FORM_DIR   FormTracer QCD project (has QCD.m)  (default below)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QCD_FORM_DIR="${QCD_FORM_DIR:-$HOME/Documents/Uni/Code/qcd-codes/vacuum/QCD_vacuum_147_AAqbq}"
OUT="$HERE/generated/form"; mkdir -p "$OUT"

if [[ "${KEEP_CACHE:-0}" != "1" ]]; then
  echo "[form] emptying ZA trace cache (FlowAA*) so the trace is recomputed ..."
  rm -f "$QCD_FORM_DIR/TraceCache/"FlowAA* 2>/dev/null || true
fi

echo "[form] regenerating in $QCD_FORM_DIR via QCD.m ..."
t0=$(date +%s.%N)
( cd "$QCD_FORM_DIR" && wolfram -script QCD.m ) 2>&1 | tee "$HERE/form_gen.log" || true
t1=$(date +%s.%N)
gen_s=$(echo "$t1 - $t0" | bc)

cp "$QCD_FORM_DIR/flows/ZA/kernel.hh" "$OUT/kernel.hh"
kb=$(( $(stat -c%s "$OUT/kernel.hh") / 1024 ))
printf "[form] gen(full QCD.m, ZA re-traced)=%.1fs  emitted=%dKB  -> generated/form/\n" "$gen_s" "$kb"
echo "GENSTAT path=form gen_s=$gen_s kb=$kb (full-suite wall; ZA MakeKernel time in form_gen.log)" >> "$HERE/bench_report.txt"
