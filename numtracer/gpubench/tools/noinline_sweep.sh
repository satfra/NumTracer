#!/usr/bin/env bash
#
# noinline_sweep.sh -- what the __noinline__ size gate does to a flow's registers and spill.
#
# Rewrites a flow's kernels.hh in place with noinline.py at each threshold, recompiles the real map
# translation unit for the target arch, and reads REG/STACK off ptxas. **Always restores the
# original header**: it is backed up before the first rewrite and restored by an EXIT trap, and the
# restore is verified with cmp.
#
# Why re-run this at all, given docs/sessions/A100_NOTES.md §5 called the gate dead: that verdict was taken when
# every kernel also carried a ~4 KB by-value interpolator copy, which dominated local traffic and
# would have hidden a spill change. With that frame gone (docs/NUMTRACER_PER_THREAD_FRAME.md),
# spill is the only local traffic left -- and on the meson codes it is 3-25x larger than on the nf2
# flows the gate was originally judged on (with_mesons ZA4: 31,832 B/thread).
#
# Usage: noinline_sweep.sh <app_dir> <flow> [thresholds...]
#   ARCH=sm_80 (default), LOGDIR=/tmp/noinline_sweep_logs
set -euo pipefail

app=${1:?app dir (the one holding flows/ and build/)}
flow=${2:?flow name, e.g. ZA4}
shift 2
thresholds=("$@")
[ ${#thresholds[@]} -eq 0 ] && thresholds=(300 500)

ARCH=${ARCH:-sm_80}
# nvcc_wrapper defaults to the system g++, which here is GCC 16 and unsupported by CUDA 12.9.
# The app build gets away with it via its own toolchain setup; a standalone replay does not.
CCBIN=${CCBIN:-/usr/bin/g++-14}
export NVCC_WRAPPER_DEFAULT_COMPILER=$CCBIN
logdir=${LOGDIR:-/tmp/noinline_sweep_logs}
mkdir -p "$logdir"

hdr="$app/flows/$flow/kernels.hh"
[ -f "$hdr" ] || { echo "no such header: $hdr" >&2; exit 1; }
tools=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
tmp=$(mktemp -d -t ntbni_XXXX)
backup="$tmp/kernels.hh.orig"
cp -p "$hdr" "$backup"

restore() {
  cp -p "$backup" "$hdr"
  if cmp -s "$backup" "$hdr"; then
    echo "[restored $hdr]"
    rm -rf "$tmp"
  else
    echo "!! RESTORE FAILED for $hdr -- original preserved at $backup" >&2
  fi
}
trap restore EXIT

# The compile command contains every shell metacharacter there is; hand it back through files
# rather than through a command substitution.
python3 "$tools/_noinline_sweep_cmd.py" "$app" "$flow" "$tmp" "$ARCH"
cwd=$(cat "$tmp/cwd.txt")

measure() {                       # $1 = label
  local log="$logdir/${flow}_$1.log"
  if ! ( cd "$cwd" && NVCC_WRAPPER_DEFAULT_COMPILER="$CCBIN" bash "$tmp/build.sh" "$tmp/$1.o" ) > "$log" 2>&1; then
    echo "  $1: COMPILE FAILED (see $log)"
    grep -m3 -iE "error" "$log" | sed 's/^/      /'
    return
  fi
  local reg stack
  reg=$(grep -oE 'Used [0-9]+ registers' "$log" | grep -oE '[0-9]+' | sort -rn | head -1)
  stack=$(grep -oE '[0-9]+ bytes stack frame' "$log" | grep -oE '^[0-9]+' | sort -rn | head -1)
  printf "  %-10s REG=%-5s STACK=%-8s\n" "$1" "${reg:-?}" "${stack:-?}"
}

echo "flow=$flow arch=$ARCH  logs=$logdir  (header restored automatically)"
echo
measure baseline
for m in "${thresholds[@]}"; do
  python3 "$tools/noinline.py" "$backup" -o "$hdr" --min "$m" >/dev/null
  measure "min$m"
done
