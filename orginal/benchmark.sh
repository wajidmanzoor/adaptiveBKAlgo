#!/usr/bin/env bash
# benchmark.sh — run every PivotBK/ReorderSib combination on all graphs
# Usage: ./benchmark.sh [data_dir]
# Reference for correctness: PivotBK+ASC

set -uo pipefail
cd "$(dirname "$0")"

BIN="./bk_algorithm"
DATA="${1:-../data}"

# macOS has no built-in 'timeout'; use gtimeout (brew coreutils) if present
if command -v gtimeout &>/dev/null; then
    RUN() { gtimeout 120 "$@"; }
elif command -v timeout &>/dev/null; then
    RUN() { timeout 120 "$@"; }
else
    RUN() { "$@"; }
fi

# ── Build ─────────────────────────────────────────────────────────────────────
echo "=== Building bk_algorithm ==="
make -s || { echo "Build failed"; exit 1; }
echo "OK"
echo ""

# ── Discover graphs ───────────────────────────────────────────────────────────
GRAPHS=()
while IFS= read -r -d '' f; do
    GRAPHS+=("$f")
done < <(find "$DATA" -maxdepth 1 -name "*.txt" -print0 | sort -z)

N=${#GRAPHS[@]}
if (( N == 0 )); then echo "No .txt files in $DATA"; exit 1; fi
echo "Graphs : $N  ($DATA)"
echo ""

# ── Helpers ───────────────────────────────────────────────────────────────────
get_count() {
    local out="$1" cnt
    cnt=$(printf '%s' "$out" | grep -oE 'Cliques Found: [0-9]+' | grep -oE '[0-9]+$') || true
    if [[ -z "$cnt" ]]; then
        cnt=$(printf '%s' "$out" | grep -oE 'cliques=[0-9]+' | grep -oE '[0-9]+$') || true
    fi
    printf '%s' "${cnt:-0}"
}

get_time() {
    local out="$1" t
    t=$(printf '%s' "$out" | grep -oE 'Time: [0-9]+\.[0-9]+' | grep -oE '[0-9.]+$') || true
    if [[ -z "$t" ]]; then
        t=$(printf '%s' "$out" | grep -oE 'time=[0-9]+\.[0-9]+' | grep -oE '[0-9.]+$') || true
    fi
    printf '%s' "${t:-0}"
}

fadd() { awk "BEGIN{printf \"%.4f\", $1 + $2}"; }

# ── Reference pass (PivotBK + ASC) ───────────────────────────────────────────
echo "=== Reference pass: PivotBK + ASC ==="
REF_COUNTS=()
REF_MS=0
for g in "${GRAPHS[@]}"; do
    out=$(RUN "$BIN" "$g" 0 1 0 2>/dev/null) || out=""
    REF_COUNTS+=("$(get_count "$out")")
    REF_MS=$(fadd "$REF_MS" "$(get_time "$out")")
done
printf "Total time: %.1f ms\n\n" "$REF_MS"

# ── Combination table ─────────────────────────────────────────────────────────
#  columns: label  algo  ord  meth
COMBOS=(
    "PivotBK      ORIG    0 0 0"
    "PivotBK      ASC     0 1 0"
    "PivotBK      DES     0 2 0"
    "---"
    "BruteForce   ORIG    1 0 0"
    "BruteForce   ASC     1 1 0"
    "BruteForce   DES     1 2 0"
    "---"
    "Backtrack    ORIG    1 0 1"
    "Backtrack    ASC     1 1 1"
    "Backtrack    DES     1 2 1"
    "---"
    "Greedy       ORIG    1 0 2"
    "Greedy       ASC     1 1 2"
    "Greedy       DES     1 2 2"
    "---"
    "Bitmask      ORIG    1 0 3"
    "Bitmask      ASC     1 1 3"
    "Bitmask      DES     1 2 3"
    "---"
    "MinHS        ORIG    1 0 4"
    "MinHS        ASC     1 1 4"
    "MinHS        DES     1 2 4"
    "---"
    "Optimized    ORIG    1 0 5"
    "Optimized    ASC     1 1 5"
    "Optimized    DES     1 2 5"
)

NC=${#COMBOS[@]}
# Count real (non-separator) combos for progress
REAL_COMBOS=0
for entry in "${COMBOS[@]}"; do
    [[ "$entry" != "---" ]] && ((REAL_COMBOS++)) || true
done

echo "=== Running $REAL_COMBOS combinations × $N graphs ==="
echo ""

WALL_START=$SECONDS

# Print table header
HDR_FMT="%-18s %-5s %10s %10s\n"
ROW_FMT="%-18s %-5s %10.1f %8s\n"
printf "$HDR_FMT" "Method" "Order" "Time(ms)" "Correct"
printf "$HDR_FMT" "------------------" "-----" "----------" "----------"

run_num=0
for entry in "${COMBOS[@]}"; do
    # Separator row
    if [[ "$entry" == "---" ]]; then
        printf "%-18s %-5s %10s %8s\n" "------------------" "-----" "----------" "----------"
        continue
    fi

    read -r label order algo ord meth <<< "$entry"
    ((run_num++)) || true

    printf "  [%2d/%d] %-12s %s ...\r" "$run_num" "$REAL_COMBOS" "$label" "$order" >&2

    ms=0
    correct=0
    for ((i=0; i<N; i++)); do
        out=$(RUN "$BIN" "${GRAPHS[$i]}" "$algo" "$ord" "$meth" 2>/dev/null) || out=""
        t=$(get_time "$out")
        ms=$(fadd "$ms" "$t")
        cnt=$(get_count "$out")
        [[ "$cnt" == "${REF_COUNTS[$i]}" ]] && ((correct++)) || true
    done

    printf "$ROW_FMT" "$label" "$order" "$ms" "$correct/$N"
done

echo ""
ELAPSED=$(( SECONDS - WALL_START ))
printf "Wall time: %dm%ds\n" $(( ELAPSED/60 )) $(( ELAPSED%60 ))
printf "Reference (PivotBK+ASC): %.1f ms — %d/%d correct by definition\n" "$REF_MS" "$N" "$N"
