#!/usr/bin/env bash
# run_methods.sh
#
# Baseline : PivotBK ASC  (clique-count reference)
# Compare  : ReorderSib × {Optimized, Backtracking} × {ORIG, ASC, DES}
#            All rules ON (prune3 always on)
#
# Usage: ./run_methods.sh [data_dir]

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"
cd "$SCRIPT_DIR"

BIN="$BK_BIN"
DATA="${1:-$DEFAULT_DATA_DIR}"

setup_timeout_runner 120

# ── Build ────────────────────────────────────────────────────────────────────
build_bk_algorithm \
  && echo "Build OK" || { echo "Build FAILED"; exit 1; }

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

# ── run_method(meth, ord) → "time_ms correct" ────────────────────────────────
# meth : 0=Backtracking  1=Optimized
# ord  : 0=ORIG  1=ASC  2=DES
run_method() {
    local meth="$1" ord="$2"
    local ms=0 correct=0
    for ((i = 0; i < N; i++)); do
        out=$(RUN "$BIN" "${GRAPHS[$i]}" 1 "$ord" "$meth" 2>/dev/null) || out=""
        ms=$(fadd "$ms" "$(get_time "$out")")
        c=$(get_cliques "$out")
        [[ "${c:-0}" == "${REF[$i]}" ]] && ((correct++)) || true
    done
    printf "%s %d" "$ms" "$correct"
}

# ── Table ────────────────────────────────────────────────────────────────────
HFMT="%-15s  %12s  %-9s  %12s  %-9s  %12s  %-9s"
RFMT="%-15s  %12s  %-9s  %12s  %-9s  %12s  %-9s"
SEP=$(printf '%-15s  %12s  %-9s  %12s  %-9s  %12s  %-9s' \
      '---------------' '------------' '---------' \
      '------------' '---------' '------------' '---------')

echo ""
printf "$HFMT\n" "Method" "ORIG ms" "ORIG ok" "ASC ms" "ASC ok" "DES ms" "DES ok"
printf '%s\n' "$SEP"

for entry in \
    "Optimized   |1" \
    "Backtracking|0"
do
    label="${entry%%|*}"
    meth="${entry#*|}"

    printf "\r  Running %-15s ORIG ...\r" "$label" >&2
    r0=$(run_method "$meth" 0)

    printf "\r  Running %-15s ASC  ...\r" "$label" >&2
    r1=$(run_method "$meth" 1)

    printf "\r  Running %-15s DES  ...\r" "$label" >&2
    r2=$(run_method "$meth" 2)

    t0="${r0% *}"; c0="${r0#* }/$N"
    t1="${r1% *}"; c1="${r1#* }/$N"
    t2="${r2% *}"; c2="${r2#* }/$N"

    printf "$RFMT\n" "$label" "${t0} ms" "$c0" "${t1} ms" "$c1" "${t2} ms" "$c2"
done

printf '%s\n' "$SEP"
printf "\nBaseline PivotBK ASC : %.3f ms   %d/%d correct by definition\n" \
       "$BASE_MS" "$N" "$N"
