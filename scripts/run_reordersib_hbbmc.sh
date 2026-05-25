#!/usr/bin/env bash
# run_reordersib_hbbmc.sh
#
# Runs:
#   1. ReorderSib (Optimized, ASC, all rules ON)
#   2. HBBMC / MCE
# on every graph in an adjacency-list dataset directory.
#
# Usage: ./run_reordersib_hbbmc.sh [data_dir] [plex]
#   data_dir  adjacency-list graph directory, defaults to repo_root/data
#   plex      HBBMC plex parameter, defaults to 0

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"
cd "$SCRIPT_DIR"

DATA_DIR="${1:-$DEFAULT_DATA_DIR}"
PLEX="${2:-0}"

BIN="$BK_BIN"
MCE="$HBBMC_BIN"

setup_timeout_runner 300

if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: '$DATA_DIR' is not a directory." >&2
    exit 1
fi

build_bk_algorithm && echo "Built bk_algorithm" || { echo "Failed to build bk_algorithm"; exit 1; }
if build_hbbmc; then
    echo "Built HBBMC/MCE"
else
    echo "Failed to build HBBMC/MCE" >&2
    if [[ "$(uname -m)" == "arm64" ]]; then
        echo "HBBMC in this tree uses x86-only SIMD intrinsics and does not build on arm64/macOS without a portability patch or x86 toolchain." >&2
    fi
    exit 1
fi
echo ""

TMP_EDGE_DIR="$(mktemp -d /private/tmp/hbbmc_edgelist.XXXXXX)"
trap 'rm -rf "$TMP_EDGE_DIR"' EXIT

python3 - "$DATA_DIR" "$TMP_EDGE_DIR" <<'PY'
import sys
from pathlib import Path

input_dir = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
files = sorted(p for p in input_dir.glob("*.txt") if p.is_file())

if not files:
    print(f"No .txt graph files found in '{input_dir}'.", file=sys.stderr)
    sys.exit(1)

for path in files:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    edges = set()
    for line in lines[1:]:
        parts = [int(x) for x in line.split()]
        if not parts:
            continue
        u = parts[0]
        for v in parts[1:]:
            if u == v:
                continue
            edges.add((u, v) if u < v else (v, u))

    out_path = output_dir / f"{path.name}.clean"
    with out_path.open("w") as f:
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")
PY

GRAPHS=()
while IFS= read -r -d '' f; do GRAPHS+=("$f"); done \
    < <(find "$DATA_DIR" -maxdepth 1 -name "*.txt" -print0 | sort -z)

N=${#GRAPHS[@]}
if [[ "$N" -eq 0 ]]; then
    echo "No graphs in $DATA_DIR" >&2
    exit 1
fi

get_hbbmc_count() {
    printf '%s\n' "$1" | grep -oE '#mc = [0-9]+' | grep -oE '[0-9]+$' | head -1
}

get_hbbmc_time() {
    printf '%s\n' "$1" | grep -oE 'runtime = [0-9]+(\.[0-9]+)? ms' | grep -oE '[0-9]+(\.[0-9]+)?' | head -1
}

get_rs_count() {
    printf '%s\n' "$1" | grep -oE 'cliques=[0-9]+' | grep -oE '[0-9]+$' | head -1
}

get_rs_time() {
    printf '%s\n' "$1" | grep -oE 'time=[0-9]+(\.[0-9]+)?' | grep -oE '[0-9]+(\.[0-9]+)?$' | head -1
}

PASS=0
FAIL=0
TIMEOUT=0

printf "Data dir : %s\n" "$DATA_DIR"
printf "Plex     : %s\n" "$PLEX"
printf "Graphs   : %d\n\n" "$N"

FMT="%-35s  %10s  %12s  %10s  %12s  %8s\n"
printf "$FMT" "Graph" "HBBMC" "HBBMC(ms)" "ReorderSib" "ReorderSib(ms)" "Status"
printf "$FMT" "-----------------------------------" "----------" "------------" "----------" "------------" "--------"

for g in "${GRAPHS[@]}"; do
    name="$(basename "$g")"
    edge_graph="$TMP_EDGE_DIR/$name.clean"

    hbbmc_out=$(RUN "$MCE" "$edge_graph" "$PLEX" 2>/dev/null) && rc_h=0 || rc_h=$?
    rs_out=$(RUN "$BIN" "$g" 1 1 1 4294967295 1 1 1 1 1 1 1 1 2>/dev/null) && rc_r=0 || rc_r=$?

    if [[ "$rc_h" -eq 124 || "$rc_r" -eq 124 ]]; then
        printf "$FMT" "$name" "TIMEOUT" "-" "TIMEOUT" "-" "SKIP"
        ((TIMEOUT++)) || true
        continue
    fi

    hbbmc_c="$(get_hbbmc_count "$hbbmc_out")"
    hbbmc_t="$(get_hbbmc_time "$hbbmc_out")"
    rs_c="$(get_rs_count "$rs_out")"
    rs_t="$(get_rs_time "$rs_out")"

    if [[ "${hbbmc_c:-X}" == "${rs_c:-Y}" ]]; then
        status="OK"
        ((PASS++)) || true
    else
        status="FAIL"
        ((FAIL++)) || true
    fi

    printf "$FMT" \
        "$name" \
        "${hbbmc_c:-?}" \
        "${hbbmc_t:-?}" \
        "${rs_c:-?}" \
        "${rs_t:-?}" \
        "$status"
done

echo ""
printf "Results: %d OK  |  %d FAIL  |  %d TIMEOUT  (out of %d graphs)\n" \
    "$PASS" "$FAIL" "$TIMEOUT" "$N"
