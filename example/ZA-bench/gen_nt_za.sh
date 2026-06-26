#!/usr/bin/env bash
# Live-regenerate the NumTracer full-QCD AAqbq ZA kernel and drop it into generated/<dest>/,
# renaming the class + support namespace so a variant can coexist with FormTracer and with the
# frozen pre-rework variant in the harness. Reports the NumTracer generation time.
#
# Usage:  ./gen_nt_za.sh [dest] [classname] [nstag]
#   dest      subdir under generated/ to write into          (default: nt)
#   classname C++ class name to rename ZA_kernel ->           (default: ZA_nt_kernel)
#   nstag     support namespace to rename za_qcd ->           (default: za_qcd_nt)
#
# Env overrides:
#   QCD_NT_DIR   NumTracer QCD project (has gen_qcd_flows.m)  (default below)
#   GLOBAL_COLLECT  if "1", pass the reworked global-collection path (NT_GLOBAL_COLLECT=1)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${1:-nt}"; CLASS="${2:-ZA_nt_kernel}"; NSTAG="${3:-za_qcd_nt}"
QCD_NT_DIR="${QCD_NT_DIR:-$HOME/Documents/Uni/Code/qcd-codes/vacuum/QCD_NumTracer}"
OUT="$HERE/generated/$DEST"; mkdir -p "$OUT"

echo "[nt] regenerating ZA in $QCD_NT_DIR (NT_FLOWS=ZA, GlobalCollect=${GLOBAL_COLLECT:-0}) ..."
# GlobalCollect: inject SetOptions[MakeNTKernel,"GlobalCollect"->True] into a temp copy of the driver
# (fold fundamental colour numerically -> group by dressing coeff -> ~FORM-compact, ~1x runtime).
driver=gen_qcd_flows.m
if [[ "${GLOBAL_COLLECT:-0}" == "1" ]]; then
  driver=gen_qcd_flows_gc.m
  sed '/SetOptions\[MakeNTKernel, Sequence @@ \$ntShimOpts\];/a SetOptions[MakeNTKernel, "GlobalCollect" -> True];' \
    "$QCD_NT_DIR/gen_qcd_flows.m" > "$QCD_NT_DIR/$driver"
fi
t0=$(date +%s.%N)
( cd "$QCD_NT_DIR" && NT_FLOWS=ZA wolfram -script "$driver" )
t1=$(date +%s.%N)
[[ "$driver" == gen_qcd_flows_gc.m ]] && rm -f "$QCD_NT_DIR/$driver"
gen_s=$(echo "$t1 - $t0" | bc)

src="$QCD_NT_DIR/flows/ZA"
# rename class (ZA_kernel) and support namespace (za_qcd) so variants don't collide in one TU.
sed -e "s/\bza_qcd\b/$NSTAG/g" -e "s/\bZA_kernel\b/$CLASS/g" "$src/kernel.hh"  > "$OUT/kernel.hh"
sed -e "s/\bza_qcd\b/$NSTAG/g"                                "$src/kernels.hh" > "$OUT/kernels.hh"

kb=$(( ($(stat -c%s "$OUT/kernel.hh") + $(stat -c%s "$OUT/kernels.hh")) / 1024 ))
ntr=$(grep -c 'inline double tr[0-9]' "$OUT/kernels.hh" || true)
printf "[nt] %s: gen=%.1fs  emitted=%dKB  trace_fns=%s  -> generated/%s/\n" "$CLASS" "$gen_s" "$kb" "$ntr" "$DEST"
echo "GENSTAT path=$DEST gen_s=$gen_s kb=$kb trace_fns=$ntr" >> "$HERE/bench_report.txt"
