#!/usr/bin/env bash
# profile_reordersib.sh
#
# Build a profiling-enabled ReorderSib binary and run it either on a curated
# set of representative graphs or on every graph in a dataset directory,
# printing the full profiling log for each individual run and saving each log
# to scripts/profile_logs/.
#
# Usage: ./scripts/profile_reordersib.sh [data_dir] [--all]
#   data_dir  graph directory, defaults to repo_root/data
#   --all     profile every graph file in the directory instead of the curated set

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"
cd "$SCRIPT_DIR"

DATA_DIR="${1:-$DEFAULT_DATA_DIR}"
MODE="${2:-}"
BIN="$BK_PROFILE_BIN"
LOG_DIR="$SCRIPT_DIR/profile_logs"

setup_timeout_runner 1800

if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: '$DATA_DIR' is not a directory." >&2
    exit 1
fi

build_bk_algorithm_profile && echo "Built profiling bk_algorithm" || {
    echo "Failed to build profiling bk_algorithm" >&2
    exit 1
}
echo ""

mkdir -p "$LOG_DIR"

PREFERRED_GRAPHS=(
    "hyves_q"
    "usa_road_sort_q"
    "randomGraph4.txt"
    "com_dblp.txt"
)

FALLBACK_GRAPHS=(
    "Deezer.txt"
    "coauthoring.txt"
    "randomGraph5.txt"
    "GSE10158_q"
    "grqc_q"
)

GRAPHS=()

select_graph_if_present() {
    local name="$1"
    local candidate="$DATA_DIR/$name"
    if [[ -f "$candidate" ]]; then
        GRAPHS+=("$candidate")
        return 0
    fi
    return 1
}

if [[ "$MODE" == "--all" ]]; then
    while IFS= read -r -d '' f; do
        GRAPHS+=("$f")
    done < <(find "$DATA_DIR" -maxdepth 1 -type f ! -name ".*" -print0 | sort -z)
else
    for name in "${PREFERRED_GRAPHS[@]}"; do
        select_graph_if_present "$name" || true
    done

    if [[ "${#GRAPHS[@]}" -lt 4 ]]; then
        for name in "${FALLBACK_GRAPHS[@]}"; do
            if [[ "${#GRAPHS[@]}" -ge 4 ]]; then
                break
            fi
            select_graph_if_present "$name" || true
        done
    fi
fi

N=${#GRAPHS[@]}
if [[ "$N" -eq 0 ]]; then
    echo "No graph files found in $DATA_DIR" >&2
    exit 1
fi

printf "Data dir : %s\n" "$DATA_DIR"
if [[ "$MODE" == "--all" ]]; then
    printf "Target set: ALL FILES\n"
else
    printf "Target set: %s\n" "${PREFERRED_GRAPHS[*]}"
fi
printf "Graphs   : %d\n" "$N"
printf "Log dir  : %s\n\n" "$LOG_DIR"

for g in "${GRAPHS[@]}"; do
    name="$(basename "$g")"
    safe_name="${name//\//_}"
    log_path="$LOG_DIR/${safe_name}.log"

    printf "==> %s\n" "$name"
    out=$(RUN "$BIN" "$g" 1 1 1 4294967295 1 1 1 1 1 1 1 1 2>&1) && rc=0 || rc=$?

    {
        printf "Graph: %s\n" "$g"
        printf "Exit code: %s\n\n" "$rc"
        printf "%s\n" "$out"
    } > "$log_path"

    printf "%s\n" "$out"
    printf "\nSaved log: %s\n\n" "$log_path"
done
