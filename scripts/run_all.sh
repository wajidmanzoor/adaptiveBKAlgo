#!/usr/bin/env bash
# run_all.sh — benchmark all method/order combinations on every graph
# Usage:  ./run_all.sh [data_dir]
# Output: printed table + results.txt

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"
cd "$SCRIPT_DIR"

BIN="$BK_BIN"
DATA="${1:-$DEFAULT_DATA_DIR}"
OUT="results.txt"

setup_timeout_runner 120

# ── Build ─────────────────────────────────────────────────────────────────────
build_bk_algorithm \
  && echo "Build OK" || { echo "Build FAILED"; exit 1; }
echo ""

# ── Graphs ────────────────────────────────────────────────────────────────────
GRAPHS=()
while IFS= read -r -d '' f; do GRAPHS+=("$f"); done \
    < <(find "$DATA" -maxdepth 1 -name "*.txt" -print0 | sort -z)
N=${#GRAPHS[@]}
[ "$N" -eq 0 ] && { echo "No graphs found in $DATA"; exit 1; }
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

# ── Reference: PivotBK + ASC ──────────────────────────────────────────────────
printf "Collecting reference (PivotBK+ASC)..."
REF_COUNTS=(); REF_MS=0
for g in "${GRAPHS[@]}"; do
    out=$(RUN "$BIN" "$g" 0 1 0 2>/dev/null) || out=""
    REF_COUNTS+=("$(get_count "$out")")
    REF_MS=$(fadd "$REF_MS" "$(get_time "$out")")
done
printf " %.1f ms\n\n" "$REF_MS"

# ── Combinations ──────────────────────────────────────────────────────────────
# Format: "Label  Order  algo ord meth"
COMBOS=(
    "PivotBK       ORIG   0 0 0"
    "PivotBK       ASC    0 1 0"
    "PivotBK       DES    0 2 0"
    "---"
    "Backtrack     ORIG   1 0 0"
    "Backtrack     ASC    1 1 0"
    "Backtrack     DES    1 2 0"
    "---"
    "Optimized     ORIG   1 0 1"
    "Optimized     ASC    1 1 1"
    "Optimized     DES    1 2 1"
)

# Count real combos for progress
TOTAL=0
for e in "${COMBOS[@]}"; do [[ "$e" != "---" ]] && ((TOTAL++)) || true; done

HDR="%-16s %-5s %10s %10s"
ROW="%-16s %-5s %10.1f %8s"
SEP="$(printf '%-16s %-5s %10s %10s' '----------------' '-----' '----------' '----------')"

# Write to both stdout and results.txt
{
    printf "$HDR\n" "Method" "Order" "Time(ms)" "Correct"
    printf '%s\n' "$SEP"
} | tee "$OUT"

run_num=0
for entry in "${COMBOS[@]}"; do
    if [[ "$entry" == "---" ]]; then
        printf '%s\n' "$SEP" | tee -a "$OUT"
        continue
    fi

    read -r label order algo ord meth <<< "$entry"
    ((run_num++)) || true
    printf "  [%2d/%d] %-12s %s ...\r" "$run_num" "$TOTAL" "$label" "$order" >&2

    ms=0; correct=0
    for ((i=0; i<N; i++)); do
        out=$(RUN "$BIN" "${GRAPHS[$i]}" "$algo" "$ord" "$meth" 2>/dev/null) || out=""
        ms=$(fadd "$ms" "$(get_time "$out")")
        [[ "$(get_count "$out")" == "${REF_COUNTS[$i]}" ]] && ((correct++)) || true
    done

    printf "$ROW\n" "$label" "$order" "$ms" "$correct/$N" | tee -a "$OUT"
done

{
    printf '%s\n' "$SEP"
    printf "\nReference (PivotBK+ASC): %.1f ms — %d/%d correct by definition\n" "$REF_MS" "$N" "$N"
    printf "Results saved to: %s/%s\n" "$(pwd)" "$OUT"
} | tee -a "$OUT"
