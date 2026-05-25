#!/usr/bin/env bash
# run_reordersib_hbbmc.sh
#
# Runs:
#   1. ReorderSib (Optimized, ASC, all rules ON)
#   2. HBBMC / MCE
# on matching graphs from two dataset directories.
#
# Usage: ./run_reordersib_hbbmc.sh <adj_dir> <edge_dir> [plex]
#   adj_dir   adjacency-list graph directory for ReorderSib
#   edge_dir  edge-list graph directory for HBBMC
#   plex      HBBMC plex parameter, defaults to 0

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"
cd "$SCRIPT_DIR"

usage() {
    echo "Usage: $0 <adj_dir> <edge_dir> [plex]"
    echo ""
    echo "  adj_dir   adjacency-list graph directory for ReorderSib"
    echo "  edge_dir  edge-list graph directory for HBBMC"
    echo "  plex      HBBMC plex parameter, default 0"
    exit 1
}

[[ $# -lt 2 ]] && usage

ADJ_DIR="$1"
EDGE_DIR="$2"
PLEX="${3:-0}"

BIN="$BK_BIN"
MCE="$HBBMC_BIN"

setup_timeout_runner 300

if [[ ! -d "$ADJ_DIR" ]]; then
    echo "Error: adjacency-list directory '$ADJ_DIR' is not a directory." >&2
    exit 1
fi

if [[ ! -d "$EDGE_DIR" ]]; then
    echo "Error: edge-list directory '$EDGE_DIR' is not a directory." >&2
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

GRAPHS=()
while IFS= read -r -d '' f; do GRAPHS+=("$f"); done \
    < <(find "$ADJ_DIR" -maxdepth 1 -type f ! -name ".*" -print0 | sort -z)

N=${#GRAPHS[@]}
if [[ "$N" -eq 0 ]]; then
    echo "No graph files found in $ADJ_DIR" >&2
    exit 1
fi

resolve_edge_graph() {
    local adj_name="$1"
    local base="${adj_name%.txt}"
    local candidates=(
        "$EDGE_DIR/$adj_name.clean"
        "$EDGE_DIR/$base.clean"
        "$EDGE_DIR/$adj_name"
        "$EDGE_DIR/$base"
    )
    local candidate
    for candidate in "${candidates[@]}"; do
        if [[ -f "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

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

format_cell() {
    local parsed="$1"
    local rc="$2"
    if [[ -n "$parsed" ]]; then
        printf '%s' "$parsed"
    else
        printf 'ERR(%s)' "$rc"
    fi
}

first_diagnostic_line() {
    printf '%s\n' "$1" | sed '/^[[:space:]]*$/d' | head -1
}

PASS=0
FAIL=0
TIMEOUT=0
MISSING=0

printf "Adj dir  : %s\n" "$ADJ_DIR"
printf "Edge dir : %s\n" "$EDGE_DIR"
printf "Plex     : %s\n" "$PLEX"
printf "Graphs   : %d\n\n" "$N"

FMT="%-35s  %10s  %12s  %10s  %12s  %8s\n"
printf "$FMT" "Graph" "HBBMC" "HBBMC(ms)" "ReorderSib" "ReorderSib(ms)" "Status"
printf "$FMT" "-----------------------------------" "----------" "------------" "----------" "------------" "--------"

for g in "${GRAPHS[@]}"; do
    name="$(basename "$g")"
    edge_graph="$(resolve_edge_graph "$name")" || edge_graph=""

    if [[ -z "$edge_graph" ]]; then
        printf "$FMT" "$name" "MISSING" "-" "-" "-" "SKIP"
        ((MISSING++)) || true
        continue
    fi

    hbbmc_out=$(RUN "$MCE" "$edge_graph" "$PLEX" 2>&1) && rc_h=0 || rc_h=$?
    rs_out=$(RUN "$BIN" "$g" 1 1 1 4294967295 1 1 1 1 1 1 1 1 2>&1) && rc_r=0 || rc_r=$?

    if [[ "$rc_h" -eq 124 || "$rc_r" -eq 124 ]]; then
        printf "$FMT" "$name" "TIMEOUT" "-" "TIMEOUT" "-" "SKIP"
        ((TIMEOUT++)) || true
        continue
    fi

    hbbmc_c="$(get_hbbmc_count "$hbbmc_out")"
    hbbmc_t="$(get_hbbmc_time "$hbbmc_out")"
    rs_c="$(get_rs_count "$rs_out")"
    rs_t="$(get_rs_time "$rs_out")"

    if [[ -z "$hbbmc_c" || -z "$hbbmc_t" || -z "$rs_c" || -z "$rs_t" ]]; then
        printf 'Parse warning for %s\n' "$name" >&2
        printf '  HBBMC rc=%s: %s\n' "$rc_h" "$(first_diagnostic_line "$hbbmc_out")" >&2
        printf '  ReorderSib rc=%s: %s\n' "$rc_r" "$(first_diagnostic_line "$rs_out")" >&2
    fi

    if [[ "${hbbmc_c:-X}" == "${rs_c:-Y}" ]]; then
        status="OK"
        ((PASS++)) || true
    else
        status="FAIL"
        ((FAIL++)) || true
    fi

    printf "$FMT" \
        "$name" \
        "$(format_cell "$hbbmc_c" "$rc_h")" \
        "$(format_cell "$hbbmc_t" "$rc_h")" \
        "$(format_cell "$rs_c" "$rc_r")" \
        "$(format_cell "$rs_t" "$rc_r")" \
        "$status"
done

echo ""
printf "Results: %d OK  |  %d FAIL  |  %d TIMEOUT  (out of %d graphs)\n" \
    "$PASS" "$FAIL" "$TIMEOUT" "$N"
if [[ "$MISSING" -gt 0 ]]; then
    printf "Skipped due to missing edge-list match: %d\n" "$MISSING"
fi
