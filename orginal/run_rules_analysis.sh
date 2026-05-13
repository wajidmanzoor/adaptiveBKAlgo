#!/usr/bin/env bash
# run_rules_analysis.sh
#
# Baseline : PivotBK ASC  (clique-count reference)
# Target   : Optimized (mode=1, meth=5) × ORIG / ASC / DES
#
# Two sections:
#   1. All 9 rules ON
#   2. One rule ON at a time (all others OFF)
#
# Usage: ./run_rules_analysis.sh [data_dir]

set -uo pipefail
cd "$(dirname "$0")"

BIN="./bk_algorithm"
DATA="${1:-../data}"

if   command -v gtimeout &>/dev/null; then RUN() { gtimeout 120 "$@"; }
elif command -v  timeout &>/dev/null; then RUN() { timeout  120 "$@"; }
else                                        RUN() {               "$@"; }
fi

# ── Build ────────────────────────────────────────────────────────────────────
make -s && echo "Build OK" || { echo "Build FAILED"; exit 1; }

# ── Graph list ───────────────────────────────────────────────────────────────
GRAPHS=()
while IFS= read -r -d '' f; do GRAPHS+=("$f"); done \
    < <(find "$DATA" -maxdepth 1 -name "*.txt" -print0 | sort -z)
N=${#GRAPHS[@]}
[ "$N" -eq 0 ] && { echo "No graphs in $DATA"; exit 1; }
echo "Graphs : $N   ($DATA)"
echo ""

# ── Helpers ──────────────────────────────────────────────────────────────────
get_cliques() {
    printf '%s\n' "$1" \
        | grep -oE 'Cliques Found: [0-9]+|cliques=[0-9]+' \
        | grep -oE '[0-9]+$' | head -1
}
get_time() {
    printf '%s\n' "$1" \
        | grep -oE 'Time: [0-9]+\.[0-9]+|time=[0-9]+\.[0-9]+' \
        | grep -oE '[0-9.]+$' | head -1
}
fadd() { awk "BEGIN{printf \"%.3f\", $1 + $2}"; }

# ── Baseline : PivotBK ASC ───────────────────────────────────────────────────
printf "Collecting baseline (PivotBK ASC) ..."
REF=()
BASE_MS=0
for g in "${GRAPHS[@]}"; do
    out=$(RUN "$BIN" "$g" 0 1 0 2>/dev/null) || out=""
    REF+=("$(get_cliques "$out")")
    BASE_MS=$(fadd "$BASE_MS" "$(get_time "$out")")
done
printf " %.3f ms\n\n" "$BASE_MS"

# ── run_config(flags, ord) → "time_ms correct" ───────────────────────────────
# flags : "p1 p2 s1 s2 s3 s4 s5 s6"  (space-separated, unquoted → word-split)
# ord   : 0=ORIG  1=ASC  2=DES
# Note: prune3 is always ON (hardcoded), not a CLI flag.
run_config() {
    local flags="$1" ord="$2"
    local ms=0 correct=0
    for ((i = 0; i < N; i++)); do
        out=$(RUN "$BIN" "${GRAPHS[$i]}" 1 "$ord" 5 4294967295 \
              $flags 2>/dev/null) || out=""
        ms=$(fadd "$ms" "$(get_time "$out")")
        c=$(get_cliques "$out")
        [[ "${c:-0}" == "${REF[$i]}" ]] && ((correct++)) || true
    done
    printf "%s %d" "$ms" "$correct"
}

# ── Config definitions ───────────────────────────────────────────────────────
# Format: "label|p1 p2 sp1 sp2 sp3 sp4 sp5 sp6"   (prune3 always ON)
CONFIGS=(
    "All ON      |1 1 1 1 1 1 1 1"
    "─────────── |───────────────"       # separator (skipped in run)
    "Only prune1 |1 0 0 0 0 0 0 0"
    "Only prune2 |0 1 0 0 0 0 0 0"
    "Only sp1    |0 0 1 0 0 0 0 0"
    "Only sp2    |0 0 0 1 0 0 0 0"
    "Only sp3    |0 0 0 0 1 0 0 0"
    "Only sp4    |0 0 0 0 0 1 0 0"
    "Only sp5    |0 0 0 0 0 0 1 0"
    "Only sp6    |0 0 0 0 0 0 0 1"
    "No rules    |0 0 0 0 0 0 0 0"
)

# Count runnable entries for progress
TOTAL_RUNS=0
for e in "${CONFIGS[@]}"; do
    [[ "${e#*|}" == *─* ]] && continue
    ((TOTAL_RUNS++)) || true
done
TOTAL_RUNS=$(( TOTAL_RUNS * 3 ))   # × 3 orders

# ── Print table ──────────────────────────────────────────────────────────────
HFMT="%-14s  %12s  %-9s  %12s  %-9s  %12s  %-9s"
RFMT="%-14s  %12s  %-9s  %12s  %-9s  %12s  %-9s"
SEP=$(printf '%-14s  %12s  %-9s  %12s  %-9s  %12s  %-9s' \
      '--------------' '------------' '---------' \
      '------------' '---------' '------------' '---------')

echo ""
printf "$HFMT\n" "Config" "ORIG ms" "ORIG ok" "ASC ms" "ASC ok" "DES ms" "DES ok"
printf '%s\n' "$SEP"

run_num=0
for entry in "${CONFIGS[@]}"; do
    label="${entry%%|*}"
    flags="${entry#*|}"

    # separator row
    if [[ "$flags" == *─* ]]; then
        printf '%s\n' "$SEP"
        continue
    fi

    printf "\r  Running %-14s ORIG ...\r" "$label" >&2
    r0=$(run_config "$flags" 0); ((run_num++)) || true

    printf "\r  Running %-14s ASC  ...\r" "$label" >&2
    r1=$(run_config "$flags" 1); ((run_num++)) || true

    printf "\r  Running %-14s DES  ...\r" "$label" >&2
    r2=$(run_config "$flags" 2); ((run_num++)) || true

    t0="${r0% *}"; c0="${r0#* }/$N"
    t1="${r1% *}"; c1="${r1#* }/$N"
    t2="${r2% *}"; c2="${r2#* }/$N"

    printf "$RFMT\n" "$label" "${t0} ms" "$c0" "${t1} ms" "$c1" "${t2} ms" "$c2"
done

printf '%s\n' "$SEP"
printf "\nBaseline PivotBK ASC : %.3f ms   %d/%d correct by definition\n" \
       "$BASE_MS" "$N" "$N"
