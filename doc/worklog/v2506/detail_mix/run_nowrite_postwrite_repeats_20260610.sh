#!/usr/bin/env bash

set -euo pipefail

ROOT="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb"
DETAIL="$ROOT/doc/worklog/v2506/detail_mix"
export OUTDIR="$DETAIL/repeat_nowrite_compute_20260610_postwrite"

mkdir -p "$OUTDIR/logs" "$OUTDIR/controlDicts"

set +u
source "$ROOT/doc/scripts/env.sh"
set -u

{
    date -Iseconds
    printf "OUTDIR=%s\n" "$OUTDIR"
    printf "dsmcFoam+=%s\n" "$(which dsmcFoam+)"
    stat -c "%y %n" "$(which dsmcFoam+)"
    stat -c "%y %n" "$ROOT/platforms/$WM_OPTIONS/lib/libdsmcFoam+.so"
    stat -c "%y %n" "$ROOT/platforms/$WM_OPTIONS/lib/libdecompose.so"
    git -C "$ROOT" rev-parse --short HEAD
    git -C "$ROOT" status --short
} > "$OUTDIR/build_state_20260610.txt"

bash "$DETAIL/run_mix_perf_nowrite_repeats_20260610.sh"
bash "$DETAIL/run_mpi8origin_nowrite_repeats_20260610.sh"
bash "$DETAIL/run_zb_mpi8origin_nowrite_repeats_20260610.sh"
python3 "$DETAIL/parse_mix_perf_nowrite_repeats_20260610.py" "$OUTDIR"
