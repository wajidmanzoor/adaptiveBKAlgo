#!/usr/bin/env bash
# run_compare.sh — PivotBK vs Optimized, all orders, per-graph results
# Usage: ./run_compare.sh [data_dir]

set -uo pipefail
cd "$(dirname "$0")"

BIN="./bk_algorithm"
DATA="${1:-../data}"

# Progress output goes directly to terminal if available
if [ -w /dev/tty ]; then
    PROGRESS_OUT="/dev/tty"
else
    PROGRESS_OUT="/dev/null"
fi

progress_bar() {
    local current="$1"
    local total="$2"
    local graph="$3"
    local width=40

    local percent=$(( current * 100 / total ))
    local filled=$(( current * width / total ))
    local empty=$(( width - filled ))

    local bar=""
    for ((i=0; i<filled; i++)); do bar+="#"; done
    for ((i=0; i<empty; i++)); do bar+="-"; done

    printf "\r[%s] %3d%%  %d/%d  %s" "$bar" "$percent" "$current" "$total" "$graph" > "$PROGRESS_OUT"

    if [ "$current" -eq "$total" ]; then
        printf "\n" > "$PROGRESS_OUT"
    fi
}

if   command -v gtimeout &>/dev/null; then RUN() { gtimeout 120 "$@"; }
elif command -v  timeout &>/dev/null; then RUN() { timeout  120 "$@"; }
else                                        RUN() {               "$@"; }
fi

make -C "$(dirname "$0")" -s && echo "Build OK" || { echo "Build FAILED"; exit 1; }
echo ""

GRAPHS=()
while IFS= read -r -d '' f; do GRAPHS+=("$f"); done \
    < <(find "$DATA" -maxdepth 1 -name "*.txt" -print0 | sort -z)

N=${#GRAPHS[@]}
[ "$N" -eq 0 ] && { echo "No graphs in $DATA"; exit 1; }

get_cliques() {
    local out="$1" c
    c=$(printf '%s\n' "$out" | grep -oE 'Cliques Found: [0-9]+' | grep -oE '[0-9]+$') || true
    [ -z "$c" ] && c=$(printf '%s\n' "$out" | grep -oE 'cliques=[0-9]+' | grep -oE '[0-9]+$') || true
    printf '%s' "${c:-0}"
}

get_time() {
    local out="$1" t
    t=$(printf '%s\n' "$out" | grep -oE 'Time: [0-9]+\.[0-9]+' | grep -oE '[0-9.]+$') || true
    [ -z "$t" ] && t=$(printf '%s\n' "$out" | grep -oE 'time=[0-9]+\.[0-9]+' | grep -oE '[0-9.]+$') || true
    printf '%s' "${t:-0}"
}

RUNS=(
    "PivotBK   ORIG   0 0 0"
    "PivotBK   ASC    0 1 0"
    "PivotBK   DES    0 2 0"
    "Optimized ORIG   1 0 5"
    "Optimized ASC    1 1 5"
    "Optimized DES    1 2 5"
)

HDR="  %-12s %-5s %12s %10s"
ROW="  %-12s %-5s %12s %10s"
SEP="  $(printf '%-12s %-5s %12s %10s' '------------' '-----' '------------' '----------')"

completed=0

for g in "${GRAPHS[@]}"; do
    graph_name="$(basename "$g")"

    progress_bar "$completed" "$N" "Running: $graph_name"

    printf "Graph: %s\n" "$graph_name"
    printf "$HDR\n" "Method" "Order" "Time(ms)" "Cliques"
    printf "%s\n" "$SEP"

    for entry in "${RUNS[@]}"; do
        read -r label order mode ord meth <<< "$entry"

        out=$(RUN "$BIN" "$g" "$mode" "$ord" "$meth" 2>/dev/null) && rc=0 || rc=$?

        if [ "$rc" -eq 124 ]; then
            t_col="     TIMEOUT"
            c_col="   TIMEOUT"
        else
            t_col=$(printf "%12.3f" "$(get_time "$out")")
            c_col=$(printf "%10s"   "$(get_cliques "$out")")
        fi

        printf "$ROW\n" "$label" "$order" "$t_col" "$c_col"
    done

    printf "\n"

    completed=$((completed + 1))
    progress_bar "$completed" "$N" "Done: $graph_name"
done