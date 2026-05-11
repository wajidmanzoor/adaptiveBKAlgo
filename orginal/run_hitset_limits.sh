#!/usr/bin/env bash
# run_hitset_limits.sh — Optimized (ORIG/ASC/DES) × hitSetLimit (1-6 + all)
# Usage: ./run_hitset_limits.sh [data_dir]

set -uo pipefail
cd "$(dirname "$0")"

BIN="./bk_algorithm"
DATA="${1:-../data}"

if   command -v gtimeout &>/dev/null; then RUN() { gtimeout 120 "$@"; }
elif command -v  timeout &>/dev/null; then RUN() { timeout  120 "$@"; }
else                                        RUN() {               "$@"; }
fi

# ── Build ─────────────────────────────────────────────────────────────────────
make -C "$(dirname "$0")" -s && echo "Build OK" || { echo "Build FAILED"; exit 1; }
echo ""

# ── Graphs ────────────────────────────────────────────────────────────────────
GRAPHS=()
while IFS= read -r -d '' f; do GRAPHS+=("$f"); done \
    < <(find "$DATA" -maxdepth 1 -name "*.txt" -print0 | sort -z)
N=${#GRAPHS[@]}
[ "$N" -eq 0 ] && { echo "No graphs in $DATA"; exit 1; }
echo "Graphs : $N  ($DATA)"
echo ""

# ── Helpers ───────────────────────────────────────────────────────────────────
get_count() {
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

fadd() { awk "BEGIN{printf \"%.4f\", $1 + $2}"; }

# ── Reference ─────────────────────────────────────────────────────────────────
printf "Collecting reference (PivotBK+ASC)..."
REF_COUNTS=(); REF_MS=0
for g in "${GRAPHS[@]}"; do
    out=$(RUN "$BIN" "$g" 0 1 0 2>/dev/null) || out=""
    REF_COUNTS+=("$(get_count "$out")")
    REF_MS=$(fadd "$REF_MS" "$(get_time "$out")")
done
printf " %.1f ms\n\n" "$REF_MS"

# ── Run: Optimized + ASC only, sweep hitSetLimit ─────────────────────────────
LIMITS=(1 2 3 4 5 6 "all")
TOTAL=${#LIMITS[@]}

printf "%-10s  %10s  %10s\n" "hitSetLimit" "Time(ms)" "Correct"
printf "%-10s  %10s  %10s\n" "----------" "----------" "----------"

for run_num in "${!LIMITS[@]}"; do
    lim=${LIMITS[$run_num]}
    printf "  [%d/%d]\r" "$((run_num+1))" "$TOTAL" >&2

    ms=0; correct=0
    for ((i=0; i<N; i++)); do
        if [[ "$lim" == "all" ]]; then
            out=$(RUN "$BIN" "${GRAPHS[$i]}" 1 1 5 2>/dev/null) || out=""
        else
            out=$(RUN "$BIN" "${GRAPHS[$i]}" 1 1 5 "$lim" 2>/dev/null) || out=""
        fi
        ms=$(fadd "$ms" "$(get_time "$out")")
        [[ "$(get_count "$out")" == "${REF_COUNTS[$i]}" ]] && ((correct++)) || true
    done

    printf "%-10s  %10.1f  %7d/%d\n" "$lim" "$ms" "$correct" "$N"
done

printf "\nReference (PivotBK+ASC): %.1f ms\n" "$REF_MS"
