#!/usr/bin/env bash
# run_per_graph.sh
#
# For each graph:
#   1. Baseline : PivotBK ASC  (clique-count + time)
#   2. ReorderSib Optimized × ORIG / ASC / DES
#      - All rules ON
#      - One rule ON at a time (others OFF)
#      - No rules
#
# Usage: ./run_per_graph.sh [data_dir]

set -uo pipefail
cd "$(dirname "$0")"

BIN="./bk_algorithm"
DATA="${1:-../data}"
SUBSET="${2:-0}"

SUBSET_NAMES=(
    "as-733"
    "as_caida"
    "AutonomousSystems"
    "ca_hepth"
    "coauthoring"
    "Deezer"
    "netscience"
    "randomGraph1"
    "randomGraph3"
    "randomGraph5"
    "yeast"
)

if   command -v gtimeout &>/dev/null; then RUN() { gtimeout 120 "$@"; }
elif command -v  timeout &>/dev/null; then RUN() { timeout  120 "$@"; }
else                                        RUN() {               "$@"; }
fi

# ── Build ────────────────────────────────────────────────────────────────────
make -s && echo "Build OK" || { echo "Build FAILED"; exit 1; }

# ── Graph list ───────────────────────────────────────────────────────────────
GRAPHS=()
if [[ "$SUBSET" == "1" ]]; then
    for name in "${SUBSET_NAMES[@]}"; do
        f="$DATA/${name}.txt"
        if [[ -f "$f" ]]; then
            GRAPHS+=("$f")
        else
            echo "WARNING: $f not found, skipping" >&2
        fi
    done
else
    while IFS= read -r -d '' f; do GRAPHS+=("$f"); done \
        < <(find "$DATA" -maxdepth 1 -name "*.txt" -print0 | sort -z)
fi
N=${#GRAPHS[@]}
[ "$N" -eq 0 ] && { echo "No graphs in $DATA"; exit 1; }
echo "Graphs : $N   ($DATA)$( [[ "$SUBSET" == "1" ]] && echo " [subset mode]" )"
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

# ── Config definitions ───────────────────────────────────────────────────────
# Format: "label|p1 p2 sp1 sp2 sp3 sp4 sp5 sp6"   (prune3 always ON)
CONFIGS=(
    "All ON      |1 1 1 1 1 1 1 1"
    "─────────── |───────────────"
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

# ── Column formats ───────────────────────────────────────────────────────────
HFMT="%-14s  %10s  %8s  %10s  %8s  %10s  %8s"
RFMT="%-14s  %10s  %8s  %10s  %8s  %10s  %8s"
SEP=$(printf '%-14s  %10s  %8s  %10s  %8s  %10s  %8s' \
      '--------------' '----------' '--------' \
      '----------' '--------' '----------' '--------')

# ── Per-graph loop ────────────────────────────────────────────────────────────
for g in "${GRAPHS[@]}"; do
    gname=$(basename "$g")

    # Baseline
    base_out=$(RUN "$BIN" "$g" 0 1 0 2>/dev/null) || base_out=""
    base_cliques=$(get_cliques "$base_out")
    base_time=$(get_time   "$base_out")

    echo "════════════════════════════════════════════════════════════════"
    echo "  Graph   : $gname"
    printf "  Baseline: PivotBK ASC → %s cliques  %s ms\n" \
           "${base_cliques:-?}" "${base_time:-?}"
    echo ""
    printf "$HFMT\n" "Config" "ORIG ms" "ORIG cliq" "ASC ms" "ASC cliq" "DES ms" "DES cliq"
    printf '%s\n' "$SEP"

    for entry in "${CONFIGS[@]}"; do
        label="${entry%%|*}"
        flags="${entry#*|}"

        if [[ "$flags" == *─* ]]; then
            printf '%s\n' "$SEP"
            continue
        fi

        # ORIG
        o0=$(RUN "$BIN" "$g" 1 0 5 4294967295 $flags 2>/dev/null) || o0=""
        t0=$(get_time "$o0"); c0=$(get_cliques "$o0")

        # ASC
        o1=$(RUN "$BIN" "$g" 1 1 5 4294967295 $flags 2>/dev/null) || o1=""
        t1=$(get_time "$o1"); c1=$(get_cliques "$o1")

        # DES
        o2=$(RUN "$BIN" "$g" 1 2 5 4294967295 $flags 2>/dev/null) || o2=""
        t2=$(get_time "$o2"); c2=$(get_cliques "$o2")

        printf "$RFMT\n" \
               "$label" \
               "${t0:-?} ms" "${c0:-?}" \
               "${t1:-?} ms" "${c1:-?}" \
               "${t2:-?} ms" "${c2:-?}"
    done

    printf '%s\n' "$SEP"
    echo ""
done

echo "Done."
